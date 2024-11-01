
    if 0:
        #the result from the programatic build
        shifter = np.array([[[-0.5, -0.5, -0.5, -0.5,  0.5,  0.5,  0.5,  0.5]],
                            [[-0.5, -0.5,  0.5,  0.5, -0.5, -0.5,  0.5,  0.5]],
                            [[-0.5,  0.5, -0.5,  0.5, -0.5,  0.5, -0.5,  0.5]]])
    if 0:
        #build
        shift_array=[-0.5, 0.5]
        shifter = None
        xyz = xyz-proj_center
        for sX in shift_array:
            for sY in shift_array:
                for sZ in shift_array:
                    this_shift = np.array([sX,sY,sZ])
                    this_shift.shape = this_shift.size,1
                    if shifter is None:
                        shifter = this_shift
                    else:
                        shifter = np.hstack([shifter, this_shift])
        shifter.shape = 3,1,8
