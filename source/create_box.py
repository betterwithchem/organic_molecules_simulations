import sys

def create_box(path_input_struct,path_box,side,shape='cubic'):

    os.system("gmx editconf"
              "-f {0}"
              "-o {1}"
              "-bt {2}"
              "-c -d 1".format(
                  path_input_struct,
                  path_box,
                  shape)
              )



    
    

    

    
