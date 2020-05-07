    t_1 = list(set(tri) - set(common)) + [common[0]] + list(set(triangles[bad_neighbor]) - set(common))
                                t_2 = list(set(triangles[bad_neighbor]) - set(common)) + [common[1]] + list(set(tri) - set(common))

                                bad_triangles[bad_neighbor] = 1
                                bad_triangles[tri_n] = 1

                                sides_1 = [tri_n] + [list(set(neighbors[bad_neighbor]) - set([tri_n]))[0]] + \
                                    [list(set(neighbors[tri_n]) - set([bad_neighbor]))[1]]
                                sides_2 = [bad_neighbor] + [list(set(neighbors[tri_n]) - set([bad_neighbor]))[0]] + \
                                    [list(set(neighbors[bad_neighbor]) - set([tri_n]))[1]]

                                for j in range(3):
                                    if neighbors[list(set(neighbors[bad_neighbor]) - set([tri_n]))[1], j] == bad_neighbor:
                                        neighbors[list(set(neighbors[bad_neighbor]) - set([tri_n]))[1], j] = tri_n
                                    if neighbors[list(set(neighbors[tri_n]) - set([bad_neighbor]))[1], j] == tri_n:
                                        neighbors[list(set(neighbors[tri_n]) - set([bad_neighbor]))[1], j] = bad_neighbor

                                neighbors[tri_n] = np.array(sides_1)
                                neighbors[bad_neighbor] = np.array(sides_2)
                                triangles[bad_neighbor] = np.array(t_1)
                                triangles[tri_n] = np.array(t_2)
                                return True