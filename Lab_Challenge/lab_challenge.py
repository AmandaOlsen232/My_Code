import math
import json
import numpy as np
import matplotlib.pyplot as plt
import MyModules as my

class CFD:
    def __init__(self, input_file):
        with open(input_file, "r") as file:
            data = json.load(file)
        #clear the output file to start
        open('mesh_NACA2412_inv.su2', 'w').close()
        
        self.plot_settings = my.MyPlot()  
        self.airfoil = data["airfoil"]["type"]
        self.open_trailing_edge = data["airfoil"]["open_trailing_edge"]
        self.r_dist = data["mesh"]["radial_distance"]
        self.wake_cells = data["mesh"]["wake_cells"]
        self.r_growth = data["mesh"]["radial_growth"]
        self.s_cells = data["mesh"]["surface_cells"]
        self.r_cells = data["mesh"]["radial_cells"]
        
        with open('mesh_NACA2412_inv.su2', 'a') as file:
            # Write the integers separated by spaces
            file.write("NDIME= 2\n")
            
    #get array describing the points of the airfoil geometry
    def geom(self):
        c = 1.0
        theta_list = np.linspace(0, 2*math.pi, self.s_cells)
        points = []
        for theta in theta_list:
            x = 0.5*(1 + math.cos(theta))
            
            ymc = int(self.airfoil[0])/100
            xmc = int(self.airfoil[1])/10
            tm = int(self.airfoil[2:])/100
            
            #determine thickness for either an open or closed trailing edge depending on user input
            if self.open_trailing_edge == True:
                t = tm*(2.969*math.sqrt(x) - 1.260*(x) - 3.516*(x)**2 + 2.843*(x)**3 - 1.015*(x)**4)
            else:
                t = tm*(2.980*math.sqrt(x) - 1.320*(x) - 3.286*(x)**2 + 2.441*(x)**3 - 0.815*(x)**4)
            
            #calculate yc and dyc_dx
            if 0 <= x <= xmc:
                yc = ymc*(2*(x/xmc) - (x/xmc)**2)
                dyc_dx = ymc*((2/xmc) - (2*x)/(xmc**2))
            elif xmc <= x <= c:
                yc = ymc*(2*((c - x)/(c - xmc)) - ((c - x)/(c - xmc))**2) 
                dyc_dx = ymc*(2*(c - x)/((c - xmc)**2) - 2/(c - xmc))
            
            #solve for either upper or lower point depending on current theta
            if theta <= math.pi:
                x_point =  x + (t/(2*math.sqrt(1 + dyc_dx**2)))*dyc_dx
                y_point = yc - (t/(2*math.sqrt(1 + dyc_dx**2)))
            elif theta >= math.pi:
                x_point = x - (t/(2*math.sqrt(1 + dyc_dx**2)))*dyc_dx
                y_point = yc + (t/(2*math.sqrt(1 + dyc_dx**2)))
            
            points.append([x_point, y_point])
            
        return np.array(points)
    
    #get distances between points
    def distances(self, growth, cells, dist):
        count = 0
        val = 0
        for _ in range(cells):
            val = val + growth**count
            count+=1
        
        #solve for the first distance
        x1 = dist/val
        distances = [x1 * growth**i for i in range(cells)]
        return np.array(distances)
    
    #generate the points for the wake plane
    def wake_plane(self):
        distances = self.distances(self.r_growth, self.wake_cells, self.r_dist)
        
        # wake_nodes_bottom, wake_nodes_top = self.wake_mesh()
        point = 1
        all_points = []
        for dist in distances:
            point = point + dist
            all_points.append([point, 0])
        

        return np.array(all_points)#, np.array(w_nodes))
    
    #calculate the tangent vectors
    def tangent(self):
        geom = self.geom()
        count = 0
        tan_list = []
        geom_list = []
        for _ in range(self.s_cells-1):
            low = geom[count]
            high = geom[count+1]
            vec = high-low
            mag = np.linalg.norm(vec)
            tan = vec/mag
            tan_list.append([tan[0], tan[1]])
            geom_list.append([geom[count][0], geom[count][1]])
            count+=1
        
        tan_list = np.array(tan_list)
        geom_list = np.array(geom_list)
        
        return geom_list, tan_list
    
    #calculate the normal vectors
    def norm(self):
        geom_list, tans_list = self.tangent()
    
        norms = []
        for row in tans_list:
            norm = [-row[1], row[0]]
            norms.append(norm)
            
        norms = np.array(norms)
        return geom_list, norms
    
    #calculate points for the cmesh
    def cmesh(self, row="na"):
        dist = self.distances(self.r_growth, self.r_cells, self.r_dist)
        norm_geom, norms = self.norm()
        
        count = 0
        points = []
        full_index = []
        index_count = 0
        for _ in range(len(norms)):
            
            if row!="na":
                count=row
                index = np.linspace(count*(self.r_cells+1), (count*(self.r_cells+1)+(self.r_cells+1))-1, self.r_cells+1)
                
            norm_vec = norms[count]
            norm_pos = norm_geom[count]
            
            total_dist = 0
            if total_dist==0:
                points.append([norm_pos[0], norm_pos[1]])
                full_index.append(index_count)
                index_count+=1
            for val in dist:
                total_dist = total_dist+val
                d = total_dist
                
                xn = norm_vec[0]*d + norm_pos[0]
                yn = norm_vec[1]*d + norm_pos[1]
                points.append([xn, yn])
                full_index.append(index_count)
                index_count+=1
            count+=1
            if row!="na":
                return np.array(points), index
        points = np.array(points)
        full_index = np.array(full_index)
        return points, full_index
    
    #calculates the points on the straight boundary
    def mesh_bounds(self):
        c = 1.0
        
        
        cmesh_bottom = self.cmesh(row=0)[0]
        cmesh_top = self.cmesh(row=-1)[0]

        xb = cmesh_bottom[-1][0]
        yb = cmesh_bottom[-1][1]
        #calculate distance between each point on the bottom
        dist = self.distances(self.r_growth, self.wake_cells, ((self.r_dist + c) - xb))
        #calculate the point locations
        bottom_points = []
        for dist_b in dist:
            point_b = [xb + dist_b, yb]
            bottom_points.append(point_b)
            xb += dist_b
        
        xt = cmesh_top[-1][0]
        yt = cmesh_top[-1][1]
        #calculate distance between each point on the top
        dist = self.distances(self.r_growth, self.wake_cells, ((self.r_dist + c) - xt))
        #calculate the point locations
        top_points = []
        for dist_t in dist:
            point_t = [xt + dist_t, yt]
            top_points.append(point_t)
            xt += dist_t
            
        return np.array(bottom_points), np.array(top_points)
    
    #calculate the normal vectors to the wake plane
    def wake_norms(self):
        #calculate the vectors between each points
        bound_bottom, bound_top = self.mesh_bounds()
        wake = self.wake_plane()
        
        count = 0
        top_vecs = []
        bottom_vecs = []
        for _ in range(len(wake)):
            b_vec = bound_bottom[count]
            t_vec = bound_top[count]
            w_vec = wake[count]
            
            top = t_vec-w_vec
            t_mag = np.linalg.norm(top)
            top_vecs.append([top[0]/t_mag, top[1]/t_mag])
            
            bottom = b_vec - w_vec
            b_mag = np.linalg.norm(bottom)
            bottom_vecs.append([bottom[0]/b_mag, bottom[1]/b_mag])
            
            
            count+=1
        
        top_vecs = np.array(top_vecs)
        bottom_vecs = np.array(bottom_vecs)
        return bottom_vecs, top_vecs
    
    #generate the mesh for the wake
    def wake_mesh(self):
        bound_bottom, bound_top = self.mesh_bounds()
        wake = self.wake_plane()
        wake_norms_b, wake_norms_t = self.wake_norms()
        
            
        count = 0
        t_points = []
        b_points = []
        for _ in range(len(wake)):
            b_point = bound_bottom[count]
            t_point = bound_top[count]
            w_point = wake[count]
            norm_b = wake_norms_b[count]
            norm_t = wake_norms_t[count]
            
            t_dist = t_point - w_point
            b_dist = b_point - w_point
            t_mag = np.linalg.norm(t_dist)
            b_mag = np.linalg.norm(b_dist)
            
            tvec_dist = self.distances(self.r_growth, self.r_cells, t_mag)
            bvec_dist = self.distances(self.r_growth, self.r_cells, b_mag)
            
            total_dist = 0
            t_points.append(w_point)
            for val in tvec_dist:
                total_dist+=val
                xn = norm_t[0]*total_dist + w_point[0]
                yn = norm_t[1]*total_dist + w_point[1]
                t_points.append([xn, yn])
                
            total_dist = 0
            for val in bvec_dist:
                total_dist+=val
                xn = norm_b[0]*total_dist + w_point[0]
                yn = norm_b[1]*total_dist + w_point[1]
                b_points.append([xn, yn])
                
            count+=1
            
        t_points = np.array(t_points)
        b_points = np.array(b_points)
        
        #get the nodes for the top and bottom meshes starting with bottom
        b_start = len(self.cmesh()[1])
        b_nodes = []
        for _ in range(len(b_points)):
            b_nodes.append(b_start)
            b_start+=1
         
        t_start = b_start
        t_nodes = []
        for _ in range(len(t_points)):
            t_nodes.append(t_start)
            t_start+=1
            
        return (b_points, np.array(b_nodes)), (t_points, np.array(t_nodes))
    
    #sort a list of 4 points in counteclockwise order
    def sort_cc(self, points):
        #sort points by x first
        sorted_points = sorted(points, key=lambda p: (p[0]))
        left_points = sorted_points[:2]
        right_points = sorted_points[2:]
        
        #sort points by y
        left_bottom, left_top = sorted(left_points, key=lambda p: p[1])
        right_bottom, right_top = sorted(right_points, key=lambda p: p[1])
        
        counterclockwise = [left_bottom.tolist(), right_bottom.tolist(), right_top.tolist(), left_top.tolist()]
        counterclockwise = np.array(counterclockwise)
    
        return counterclockwise
    
    def add_cmesh(self):
        row_count = 0
        final_ints = []
        cell_count = 0
        #the first for loop tells the row you are on
        for _ in range(self.s_cells-2):
            row1, index1 = self.cmesh(row=row_count)
            row2, index2 = self.cmesh(row=row_count+1)
            #now I need a for loop for the point you are on
            point_count=0
            for _ in range(len(index1)-1):
                cell_count+=1
                first = row1[point_count].tolist()
                firsti = index1[point_count]
                first.append(firsti)
                
                second = row1[point_count+1].tolist()
                secondi = index1[point_count+1]
                second.append(secondi)
                
                third = row2[point_count].tolist()
                thirdi = index2[point_count]
                third.append(thirdi)
                
                fourth = row2[point_count+1].tolist()
                fourthi = index2[point_count+1]
                fourth.append(fourthi)
                
                points = [np.array(first), np.array(second), np.array(third), np.array(fourth)]
                sorted_points = self.sort_cc(points)
                plt.plot(sorted_points[:,0], sorted_points[:,1])
                ind_1 = int(sorted_points[0][2])
                ind_2 = int(sorted_points[1][2])
                ind_3 = int(sorted_points[2][2])
                ind_4 = int(sorted_points[3][2])
                
                final_ints.append([ind_1, ind_2, ind_3, ind_4])
        
                point_count+=1
            
            row_count+=1
        
        row = 0
        for _ in range(self.wake_cells-1):
            left = self.w_mesh_row_nodes(row)[0]
            right = self.w_mesh_row_nodes(row+1)[0]
            point = 0
            for _ in range(len(left[1])-1):
                cell_count+=1
                one = left[0][point].tolist()
                one.append(left[1][point])
                
                two = left[0][point+1].tolist()
                two.append(left[1][point+1])
                
                three = right[0][point].tolist()
                three.append(right[1][point])
                
                four = right[0][point+1].tolist()
                four.append(right[1][point+1])
                
                points = [np.array(one), np.array(two), np.array(three), np.array(four)]
                sorted_points = self.sort_cc(points)
                # plt.plot(sorted_points[:,0], sorted_points[:,1])
                
                ind_1 = int(sorted_points[0][2])
                ind_2 = int(sorted_points[1][2])
                ind_3 = int(sorted_points[2][2])
                ind_4 = int(sorted_points[3][2])
                
                final_ints.append([ind_1, ind_2, ind_3, ind_4])
                    
                point+=1
            row+=1
            
        #now the top mesh
        row = 0
        for _ in range(self.wake_cells-1):
            left = self.w_mesh_row_nodes(row)[1]
            right = self.w_mesh_row_nodes(row+1)[1]
            point = 0
            for _ in range(len(left[1])-1):
                cell_count+=1
                one = left[0][point].tolist()
                one.append(left[1][point])
                
                two = left[0][point+1].tolist()
                two.append(left[1][point+1])
                
                three = right[0][point].tolist()
                three.append(right[1][point])
                
                four = right[0][point+1].tolist()
                four.append(right[1][point+1])
                
                points = [np.array(one), np.array(two), np.array(three), np.array(four)]
                sorted_points = self.sort_cc(points)
                # plt.plot(sorted_points[:,0], sorted_points[:,1])
                
                ind_1 = int(sorted_points[0][2])
                ind_2 = int(sorted_points[1][2])
                ind_3 = int(sorted_points[2][2])
                ind_4 = int(sorted_points[3][2])
                
                final_ints.append([ind_1, ind_2, ind_3, ind_4])
                    
                point+=1
            row+=1
            
        row = 0
        #connect the top and bottom wake mesh
        for _ in range(self.wake_cells-1):
            cell_count+=1
            left_b, left_t = self.w_mesh_row_nodes(row)
            right_b, right_t = self.w_mesh_row_nodes(row+1)
            
            one = left_b[0][0].tolist()
            one.append(left_b[1][0])
            
            two = left_t[0][0].tolist()
            two.append(left_t[1][0])
            
            three = right_b[0][0].tolist()
            three.append(right_b[1][0])
            
            four = right_t[0][0].tolist()
            four.append(right_t[1][0])
            
            points = [np.array(one), np.array(two), np.array(three), np.array(four)]
            sorted_points = self.sort_cc(points)
            
            # plt.plot(sorted_points[:,0], sorted_points[:,1])
            ind_1 = int(sorted_points[0][2])
            ind_2 = int(sorted_points[1][2])
            ind_3 = int(sorted_points[2][2])
            ind_4 = int(sorted_points[3][2])
            
            final_ints.append([ind_1, ind_2, ind_3, ind_4])
                
            row+=1
        
        top_c = self.cmesh(row=self.s_cells-2)
        w_top = self.w_mesh_row_nodes(0)[1] 
        
        point = 0
        for _ in range(len(top_c[1])-1):
            cell_count+=1
            one = w_top[0][point].tolist()
            one.append(w_top[1][point])
            
            two = w_top[0][point+1].tolist()
            two.append(w_top[1][point+1])
            
            three = top_c[0][point].tolist()
            three.append(top_c[1][point])
            
            four = top_c[0][point+1].tolist()
            four.append(top_c[1][point+1])
            
            points = [np.array(one), np.array(two), np.array(three), np.array(four)]
            sorted_points = self.sort_cc(points)
            print(sorted_points)
            # plt.plot(sorted_points[:,0], sorted_points[:,1])
            ind_1 = int(sorted_points[0][2])
            ind_2 = int(sorted_points[1][2])
            ind_3 = int(sorted_points[2][2])
            ind_4 = int(sorted_points[3][2])
            
            final_ints.append([ind_1, ind_2, ind_3, ind_4])

            point+=1
            
        b_c = self.cmesh(row=0)
        w_b, w_t = self.w_mesh_row_nodes(0)
        vals = np.vstack((w_t[0][0], w_b[0]))
        nodes = np.vstack((np.array(w_t[1][0]).reshape(1,1), w_b[1].reshape(self.r_cells,1)))
        w_b = (vals, nodes.reshape(self.r_cells+1))
        
        point = 0
        for _ in range(len(b_c[1])-1):
            cell_count+=1
            one = w_b[0][point].tolist()
            one.append(w_b[1][point])
            
            two = w_b[0][point+1].tolist()
            two.append(w_b[1][point+1])
            
            three = b_c[0][point].tolist()
            three.append(b_c[1][point])
            
            four = b_c[0][point+1].tolist()
            four.append(b_c[1][point+1])
            
            points = [np.array(one), np.array(two), np.array(three), np.array(four)]
            sorted_points = self.sort_cc(points)
            plt.plot(sorted_points[:,0], sorted_points[:,1])
            ind_1 = int(sorted_points[0][2])
            ind_2 = int(sorted_points[1][2])
            ind_3 = int(sorted_points[2][2])
            ind_4 = int(sorted_points[3][2])
            
            final_ints.append([ind_1, ind_2, ind_3, ind_4])
            point+=1
        with open('mesh_NACA2412_inv.su2', 'a') as file:
            # Write the integers separated by spaces
            file.write(f"NELEM= {cell_count}\n")
            for row in final_ints:
                file.write(f"9 {row[0]} {row[1]} {row[2]} {row[3]}\n")

    def add_points(self):
        c = self.cmesh()
        bottom, top = self.wake_mesh()
        
        one, one_n = c
        two, two_n = bottom
        three, three_n = top
        
        one = np.hstack((one, one_n.reshape(len(one_n),1)))
        two = np.hstack((two, two_n.reshape(len(two_n),1)))
        three = np.hstack((three, three_n.reshape(len(three_n),1)))
        
        first = np.vstack((one, two))
        all_points = np.vstack((first, three))
        total_num = int(all_points[-1][2]+1)
        
        with open('mesh_NACA2412_inv.su2', 'a') as file:
            # Write the integers separated by spaces
            file.write(f"NPOIN= {total_num}\n")
            for row in all_points:
                file.write(f"{row[0]} {row[1]} {int(row[2])}\n")

    def add_sbounds(self):
        #add the airfoil boundaries
        count=0
        all_points = []
        all_nodes = []
        for _ in range(self.s_cells-1):
            s_point, s_node = self.cmesh(row=count)
            s_point = s_point[0]
            s_node = int(s_node[0])
            all_points.append([s_point[0], s_point[1]])
            all_nodes.append(s_node)
            count+=1
        all_points = np.array(all_points)
        all_nodes = np.array(all_nodes)
        
        with open('mesh_NACA2412_inv.su2', 'a') as file:
            # Write the integers separated by spaces
            file.write("NMARK= 2\n")
            file.write("MARKER_TAG= airfoil\n")
            file.write(f"MARKER_ELEMS= {len(all_nodes)}\n")
            
            count = 0
            for _ in range(len(all_nodes)-1):
                first = all_nodes[count]
                second = all_nodes[count+1]
                file.write(f"3 {first} {second}\n")
                count+=1
            
            file.write(f"3 {all_nodes[-1]} {all_nodes[0]}\n")
            
    #input a row and it tells me the nodes on that row for the wake mesh
    def w_mesh_row_nodes(self, row):
        wake_nodes_bottom, wake_nodes_top = self.wake_mesh()
        
        num_wn_t = self.r_cells + 1
        num_wn_b = self.r_cells
         
        point = 0
        row_nodes_t = []
        points_t = []
        for _ in range(num_wn_t):
            node = row*num_wn_t + point
            row_nodes_t.append(wake_nodes_top[1][node])
            points_t.append(wake_nodes_top[0][node])
            point+=1
            
        point = 0
        row_nodes_b = []
        points_b = []
        for _ in range(num_wn_b):
            node = row*num_wn_b + point
            row_nodes_b.append(wake_nodes_bottom[1][node])
            points_b.append(wake_nodes_bottom[0][node])
            point+=1
            
        return (np.array(points_b), np.array(row_nodes_b)), (np.array(points_t), np.array(row_nodes_t))
    
    #get the points for the outside boundary
    def get_outside_points(self):
        #now add outside bounds starting with cmesh
        count=0
        all_points = []
        all_nodes = []
        for _ in range(self.s_cells-1):
            s_point, s_node = self.cmesh(row=count)
            s_point = s_point[-1]
            s_node = int(s_node[-1])
            all_points.append([s_point[0], s_point[1]])
            all_nodes.append(s_node)
            count+=1
            
        #now add top row from w_mesh
        row = 0
        for _ in range(self.wake_cells):
            bottom, top = self.w_mesh_row_nodes(row)
            t_point = top[0][-1]
            t_node = top[1][-1]
            all_points.append([t_point[0], t_point[1]])
            all_nodes.append(t_node)
            row+=1
        
        #now add right top boundary
        bottom, top = self.w_mesh_row_nodes(-1)
        flipped_top = top[0][::-1]
        flipped_index = top[1][::-1]
        count=1
        for _ in range(len(flipped_top)-1):
            rt_point = flipped_top[count]
            rt_node = flipped_index[count]
            all_points.append([rt_point[0], rt_point[1]])
            all_nodes.append(rt_node)
            count+=1
            
        #now add right bottom boundary
        bottom, top = self.w_mesh_row_nodes(-1)
        count=0
        for _ in range(len(bottom[1])):
            rb_point = bottom[0][count]
            rb_node = bottom[1][count]
            all_points.append([rb_point[0], rb_point[1]])
            all_nodes.append(rb_node)
            count+=1
        
        #now add final line
        row=-2
        for _ in range(self.wake_cells-1):
            bottom, top = self.w_mesh_row_nodes(row)
            b_point = bottom[0][-1]
            b_node = bottom[1][-1]
            all_points.append([b_point[0], b_point[1]])
            all_nodes.append(b_node)
            row-=1
        
        all_points = np.array(all_points)
        all_nodes = np.array(all_nodes)
        return all_points, all_nodes
    
    #add the outside boundary data
    def add_outside_bounds(self):
        points, nodes = hello.get_outside_points()
        with open('mesh_NACA2412_inv.su2', 'a') as file:
            # Write the integers separated by spaces
            file.write("MARKER_TAG= farfield\n")
            file.write(f"MARKER_ELEMS= {len(nodes)}\n")
            
            count = 0
            for _ in range(len(points)-1):
                first = nodes[count]
                second = nodes[count+1]
                file.write(f"3 {first} {second}\n")
                count+=1
            
            file.write(f"3 {nodes[-1]} {nodes[0]}\n")

    #plot everything
    def plot(self, plot_geom=True, plot_wake=False, plot_tan=False, plot_norm=False, plot_cmesh=False, plot_wmesh=False):
        
        #plot airfoil geometry and wake plane
        if plot_geom==True:
            airfoil_geom = self.geom()
            plt.plot(airfoil_geom[:, 0], airfoil_geom[:, 1])
            
        if plot_wake==True:
            wake_plane_points = self.wake_plane()
            plt.scatter(wake_plane_points[:, 0], wake_plane_points[:, 1], s=2, color="teal")
        
        if plot_tan==True:
            geom_list, tan_list = self.tangent()
            plt.quiver(geom_list[:,0], geom_list[:,1], tan_list[:,0], tan_list[:,1], color="red")
        
        if plot_norm==True:
            geom_list, norms = self.norm()
            plt.quiver(geom_list[:,0], geom_list[:,1], norms[:,0], norms[:,1], color="blue")
        
        if plot_cmesh==True:
            cmesh_points, index = self.cmesh()
            plt.scatter(cmesh_points[:, 0], cmesh_points[:, 1], s=2, color="teal")
        
        if plot_wmesh==True:
            wmesh_b, wmesh_t = self.wake_mesh()
            plt.scatter(wmesh_b[0][:,0], wmesh_b[0][:,1], s=2, color="teal")
            plt.scatter(wmesh_t[0][:,0], wmesh_t[0][:,1], s=2, color="teal")

        plt.xlabel("x", fontstyle='italic')
        plt.ylabel("y", fontstyle='italic')
        plt.axis("equal")
        
        plt.show()
        # plt.close()
        
        
if __name__=="__main__":
    hello = CFD("input_file.json")
    hello.plot(plot_cmesh=True, plot_wmesh=True)

    hello.add_cmesh()
    hello.add_points()
    hello.add_sbounds()
    hello.add_outside_bounds()

