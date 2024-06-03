import numpy as np

##The functions below write two different .lammpstrj files

def analysis_Initial(filename,no_of_atom,boxL,InitialCentralPosition, InitialTag2pos, InitialTag3pos): 
    #Create LAMMPS output file where tag and particles are in the same file
    ##includes molecule number 
    f= open(filename,"w")
    
    
    
    f.write('ITEM: TIMESTEP')
    f.write('\n')
    f.write('%i ' %0)
    f.write('\n')
    f.write('ITEM: NUMBER OF ATOMS')
    f.write('\n')
    f.write('%i ' %(no_of_atom*3))
            
    f.write('\n')
    f.write('ITEM: BOX BOUNDS pp pp pp')
    f.write('\n')
    f.write('0  %.0f' %boxL)
    f.write('\n')
    f.write('0  %.0f' %boxL)
    f.write('\n')
    f.write('0  %.0f' %boxL)
    f.write('\n')
    f.write('ITEM: ATOMS molID id type x y z')
    f.write('\n')
            
            
            
    for i in range(no_of_atom):
        
                
            f.write(str(2*i+(i+1))+ ' %i 1 ' %(i+1)) #%i 1 gives molecule number and type
                  
            f.write('%.10f ' %InitialCentralPosition[i,0])
            f.write('%.10f ' %InitialCentralPosition[i,1])
            f.write('%.1f ' %0)                
            f.write('\n')  
    
                
            f.write(str(2*i+(i+2))+ ' %i 2 ' %(i+1))
                  
            f.write('%.10f ' %InitialTag2pos[i,0])
            f.write('%.10f ' %InitialTag2pos[i,1])
            f.write('%.10f ' %0)                
            f.write('\n')

            f.write(str(2*i+(i+3))+ ' %i 3 ' %(i+1))
                  
            f.write('%.10f ' %InitialTag3pos[i,0])
            f.write('%.10f ' %InitialTag3pos[i,1])
            f.write('%.10f ' %0)                
            f.write('\n')
                

        
    f.close()

def analysis_step(filename,no_of_atom,boxL,centralPosition,tag2pos, tag3pos,step): 
    #Append to the file created in analysis_initial()
    f = open(filename,"a")
    
    
    
    f.write('ITEM: TIMESTEP')
    f.write('\n')
    f.write('%i ' %(step+1))
    f.write('\n')
    f.write('ITEM: NUMBER OF ATOMS')
    f.write('\n')
    f.write('%i ' %(no_of_atom*3))
            
    f.write('\n')
    f.write('ITEM: BOX BOUNDS pp pp pp')
    f.write('\n')
    f.write('0  %.0f' %boxL)
    f.write('\n')
    f.write('0  %.0f' %boxL)
    f.write('\n')
    f.write('0  %.0f' %boxL)
    f.write('\n')
    f.write('ITEM: ATOMS molID id type x y z')
    f.write('\n')
            
            
            
    for i in range(no_of_atom):
        
    
        
                
            f.write(str(2*i+(i+1))+ ' %i 1 ' %(i+1)) 
                  
            f.write('%.10f ' %centralPosition[i,0])
            f.write('%.10f ' %centralPosition[i,1])
            f.write('%.1f ' %0)                
            f.write('\n')  
    
                
            f.write(str(2*i+(i+2))+ ' %i 2 ' %(i+1))
                  
            f.write('%.10f ' %tag2pos[i,0])
            f.write('%.10f ' %tag2pos[i,1])
            f.write('%.10f ' %0)                
            f.write('\n')

            f.write(str(2*i+(i+3))+ ' %i 3 ' %(i+1))
                  
            f.write('%.10f ' %tag3pos[i,0])
            f.write('%.10f ' %tag3pos[i,1])
            f.write('%.10f ' %0)                
            f.write('\n')
                

        
    f.close()

import numpy as np

def LAMMPS_Initial(filename,no_of_atom,boxL,InitialCentralPosition, InitialTag2pos, InitialTag3pos): 
    #.lammpstrj without molecule 
    f= open(filename,"w")
    
    
    
    f.write('ITEM: TIMESTEP')
    f.write('\n')
    f.write('%i ' %0)
    f.write('\n')
    f.write('ITEM: NUMBER OF ATOMS')
    f.write('\n')
    f.write('%i ' %(no_of_atom*3))
            
    f.write('\n')
    f.write('ITEM: BOX BOUNDS pp pp pp')
    f.write('\n')
    f.write('0  %.0f' %boxL)
    f.write('\n')
    f.write('0  %.0f' %boxL)
    f.write('\n')
    f.write('0  %.0f' %boxL)
    f.write('\n')
    f.write('ITEM: ATOMS id type x y z')
    f.write('\n')
            
            
            
    for i in range(no_of_atom):
        
                
            f.write(str(2*i+(i+1))+ ' 1 ' ) #%i 1 gives molecule number and type
                  
            f.write('%.10f ' %InitialCentralPosition[i,0])
            f.write('%.10f ' %InitialCentralPosition[i,1])
            f.write('%.1f ' %0)                
            f.write('\n')  
    
                
            f.write(str(2*i+(i+2))+ ' 2 ' )
                  
            f.write('%.10f ' %InitialTag2pos[i,0])
            f.write('%.10f ' %InitialTag2pos[i,1])
            f.write('%.10f ' %0)                
            f.write('\n')

            f.write(str(2*i+(i+3))+ ' 3 ' )
                  
            f.write('%.10f ' %InitialTag3pos[i,0])
            f.write('%.10f ' %InitialTag3pos[i,1])
            f.write('%.10f ' %0)                
            f.write('\n')
                

        
    f.close()

def LAMMPS_step(filename,no_of_atom,boxL,centralPosition,tag2pos, tag3pos,step): #To Create LAMMPS output filew which tag + particles are in the same file
    #.lammpstrj without molecule 
    f = open(filename,"a")
    
    
    
    f.write('ITEM: TIMESTEP')
    f.write('\n')
    f.write('%i ' %(step+1))
    f.write('\n')
    f.write('ITEM: NUMBER OF ATOMS')
    f.write('\n')
    f.write('%i ' %(no_of_atom*3))
            
    f.write('\n')
    f.write('ITEM: BOX BOUNDS pp pp pp')
    f.write('\n')
    f.write('0  %.0f' %boxL)
    f.write('\n')
    f.write('0  %.0f' %boxL)
    f.write('\n')
    f.write('0  %.0f' %boxL)
    f.write('\n')
    f.write('ITEM: ATOMS id type x y z')
    f.write('\n')
            
            
            
    for i in range(no_of_atom):
        
    
        
                
            f.write(str(2*i+(i+1))+ ' 1 ' ) 
                  
            f.write('%.10f ' %centralPosition[i,0])
            f.write('%.10f ' %centralPosition[i,1])
            f.write('%.1f ' %0)                
            f.write('\n')  
    
                
            f.write(str(2*i+(i+2))+ ' 2 ' )
                  
            f.write('%.10f ' %tag2pos[i,0])
            f.write('%.10f ' %tag2pos[i,1])
            f.write('%.10f ' %0)                
            f.write('\n')

            f.write(str(2*i+(i+3))+ ' 3 ' )
                  
            f.write('%.10f ' %tag3pos[i,0])
            f.write('%.10f ' %tag3pos[i,1])
            f.write('%.10f ' %0)                
            f.write('\n')
                

        
    f.close()


