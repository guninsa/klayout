
# Enter your Python code here

import numpy as np
from math import sqrt, cos, sin, atan2, pi, copysign

from importlib import reload
import pya
from pya import Point,DPoint,DSimplePolygon,SimplePolygon, DPolygon, Polygon,  Region
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from classLib import *
#import classLib
#reload(classLib)
#from classLib.chipHole import CirclePolygon, ScrewHoles, MetalRounds, SampleHoles
#from classLib.coplanars import CPWParameters, CPW, CPWOpened, CPWArc, CPW2CPW, CPW2CPWArc, Coil_type_1, CPWRLPath, DPathCPW, Bridge1, VIA, CirclePolygon, Holes,  Circles, CoaxTl
#reload(classLib.coplanars)
#reload(classLib.chipHole)

#mport classLib
#reload(classLib)
#from classLib.coplanars import VIA, CPW, CPW
#############

if __name__ == "__main__":
# getting main references of the application
    app = pya.Application.instance()
    mw = app.main_window()
    lv = mw.current_view()
    cv = None
    
    #this insures that lv and cv are valid objects
    if( lv == None ):
        cv = mw.create_layout(1)
        lv = mw.current_view()
    else:
        cv = lv.active_cellview()

# find or create the desired by programmer cell and layer
    layout = cv.layout()
    layout.dbu = 0.001

    cell = layout.create_cell( "testScript" )
    cell.clear()

    infos = [pya.LayerInfo(i+1,0) for i in range(6)] # range 6
    layers = [layout.layer(infos[i]) for i in range(6)] #range 6
   
   


### DRAW SECTION START ###

    def draw_chip():
      pass
      
    radius = 55e6
    nr_points = 32
    angles = np.linspace(0,2*np.pi,nr_points+1)[0:-1]
    points = []
    
    for ind,angle in enumerate(angles):
      points.append(pya.Point(radius*np.cos(angle),radius*np.sin(angle)))
    circle = pya.SimplePolygon(points)
    
    for layer in layers[0:5]:
      cell.shapes(layer).insert(circle)

  
    # setting layout view  
    lv.select_cell(cell.cell_index(), 0)
    lv.add_missing_layers()

   
    ### DRAW SECTION END ###
    
    lv.zoom_fit()
 ### DRAW SECTION START ###
    

    
    #width = 10e4
   # gap = 20e4


   # feedline = CPW(width, gap, start=DPoint(10e5,10e5), end=DPoint(200e5,200e3))
    #feedline.place(cell, layers[0]) 
    
    
    hole = ScrewHoles(100e5, origin = DPoint(300e5, 300e5))
    #hole.place(cell, layers[9])
    
    cpw_params = {'w':150e3, 's':200e3, 'g':150e3}
    sample_hole = SampleHoles(dx = 14e6, dy = 14e6, step_forward = 100e3, origin = DPoint(0,0), gap_from_step = 1e3,
               trans_in = None, ports_num = 40, ports = [], cpw_params = cpw_params, rounds_radius = 1e3, via_r = 0,
                 region_id = "default")

    
    sample_hole_stepped = SampleHoles(dx = 14e6, dy = 14e6, step_forward = 300e3, origin = DPoint(0,0), gap_from_step = -1e3,
               trans_in = None, ports_num = 40, ports = [], cpw_params = cpw_params, rounds_radius = 1e3, via_r = 100e3,
                 region_id = "default")
                

    i = 0
    for port in sample_hole_stepped.ports:
      i+= 1
      if 1<=i<=10:
          opened_test = CPWOpened(cpw_params['s'], cpw_params['g'], start=port, end=port, radius=None, number_of_points=32., trans_in = DTrans.R180 )
          #opened_test.place(cell, layers[3])
      elif 11<=i<=20:
          opened_test = CPWOpened(cpw_params['s'], cpw_params['g'], start=port, end=port, radius=None, number_of_points=32., trans_in = DTrans.R90 )
         # opened_test.place(cell, layers[3])
      elif 21<=i<=30:
          opened_test = CPWOpened(cpw_params['s'], cpw_params['g'], start=port, end=port, radius=None, number_of_points=32., trans_in = DTrans.R0 )
         # opened_test.place(cell, layers[3])
      elif 31<=i<=40:
          opened_test = CPWOpened(cpw_params['s'], cpw_params['g'], start=port, end=port, radius=None, number_of_points=32., trans_in = DTrans.R270 )
          #opened_test.place(cell, layers[3])
          
          
     # x1 = 14.6e3/2
     #y1 = -14.6e3/2
     #x2 = -14.6e3/2
     #y2 =  14.6e3/2
    
    #k = (y1-y2)/(x1-x2)
    #b = y1 - x1*(y1-y2)/(x1-x2)
    

    offset_y = 40e6
    offset_x = 25e6
    offset_step = 2*(radius - offset_x)/11
    i =0
    coaxials = []
    for port in sample_hole_stepped.ports[0:10]:
        i+=1
        coax = CoaxTl(center = DPoint(-radius + offset_x+offset_step*i,  offset_y), trans_in = DTrans.R180, cell= cell, layer1 = layers[3],  layer2 = layers[2], layer3 = layers[1],  layer4 =layers[0], via_layer = layers[5], orientation = np.pi)
        coaxials.append(coax)
        
    i =0
    for port in sample_hole_stepped.ports[10:20]:
         i+=1
         coax = CoaxTl(center = DPoint( offset_y , radius - offset_x-offset_step*i), trans_in = DTrans.R90, cell= cell, layer1 = layers[3],  layer2 = layers[2], layer3 = layers[1],  layer4 =layers[0], via_layer = layers[5], orientation = 0.5*np.pi)
         coaxials.append(coax)
         
    i =0
    for port in sample_hole_stepped.ports[30:40]:
         i+=1
         coax = CoaxTl(center = DPoint(-radius + offset_x+offset_step*i,  -offset_y), trans_in = DTrans.R0, cell= cell, layer1 = layers[3],  layer2 = layers[2], layer3 = layers[1],  layer4 =layers[0], via_layer = layers[5], orientation = 0*np.pi)
         coaxials.append(coax)

    i =0
    for port in sample_hole_stepped.ports[10:20]:
         i+=1
         coax = CoaxTl(center = DPoint( -offset_y , radius - offset_x-offset_step*i), trans_in = DTrans.R270, cell= cell, layer1 = layers[3],  layer2 = layers[2], layer3 = layers[1],  layer4 =layers[0], via_layer = layers[5], orientation = 1.5*np.pi)  
         coaxials.append(coax)
         

    cpw_params = CPWParameters(200e3, 150e3)
    
         
    for i in range(10):
  
        if i<=4:
            flag_len = 1
            flag_ang = 1
        elif i==5:
            flag_len = 0
            flag_ang = -1
        elif i>5:
            flag_len = -1
            flag_ang = -1
            
        step += flag_len*2e6

        segment_lenghts = [abs(coaxials[i].cpw_end.y-sample_hole_stepped.ports[i].y)-step-abs(coaxials[i].cpw_end.x-sample_hole_stepped.ports[i].x), abs(coaxials[i].cpw_end.x-sample_hole_stepped.ports[i].x)*np.sqrt(2),step ]
    
        feedline1 = CPWRLPath(coaxials[i].cpw_end,"LRLRL", cpw_params, 100e4, segment_lenghts, [flag_ang*pi/4,-pi/4*flag_ang], trans_in = DTrans.R270)
        feedline1.place(cell, layers[3])


        feedline2 = CPWRLPath(coaxials[i+10].cpw_end,"LRLRL", cpw_params, 100e4, segment_lenghts, [flag_ang*pi/4,-pi/4*flag_ang], trans_in = DTrans.R180)
        feedline2.place(cell, layers[3])
        
        
        feedline3 = CPWRLPath(sample_hole_stepped.ports[i+20],"LRLRL", cpw_params, 100e4, segment_lenghts[::-1], [flag_ang*pi/4,-pi/4*flag_ang], trans_in = DTrans.R270)
        feedline3.place(cell, layers[3])
        

        feedline4 = CPWRLPath(sample_hole_stepped.ports[i+30],"LRLRL", cpw_params, 100e4, segment_lenghts[::-1], [flag_ang*pi/4,-pi/4*flag_ang], trans_in = DTrans.R180)
        feedline4.place(cell, layers[3])



    

    

    
    
    
    #hole_end = ScrewHoles(10e5, origin = feedline4.end)
    #hole_end.place(cell, layers[9])
    
    #for name, primitive in feedline4.primitives.items():
      #print(name, primitive)
      
      
    
    #opened_test.place(cell, layers[1])
    
    #opened_testopened_test = CPW(11e3, 5e3, start=DPoint(0,-500e3), end=opened_test.start, trans_in = DTrans.R0)
    #opened_testopened_test.place(cell, layers[1]) 
    
    #print(np.linspace(0, np.pi, nr_points, endpoint=True))
    #print('x= ' ,DPoint(5e3,7e3).x)
    #print('x= ' ,DPoint(1e3,3e3).y)
    #print(DPoint(1e3,2e3)+DPoint(3e3,4e3))

   # via = VIA(center=DPoint(1e5,1e5), radius = 20e5, number_of_points=32, via_width = 10e5, trans_in = None)
    
    #vias = via.viafy_CPW(opened_testopened_test, 50e3,80e3,10e3, dest=cell, gnd2gnd_dy=None,  bridge_layer1 = layers[4], bridge_layer2= layers[4], dest2=None, via_left = False, via_right=False,
      #number_of_points=321, avoid_points = [], avoid_distances = [])
    #print(opened_testopened_test.dr.y/opened_testopened_test.dr.x)
    #print(isinstance(feedline4.primitives['cpw_0'], CPW))
    
    #vias2 = via.viafy_CPW(feedline4, 50e3,80e3,10e3, dest=cell, gnd2gnd_dy=None,  bridge_layer1 = layers[4], bridge_layer2= layers[4], dest2=None, via_left = False, via_right=False,
      #number_of_points=321, avoid_points = [], avoid_distances = [], flag = 2)
    #print('S' in 'S1')

#    for name, primitive in feedline4.primitives.items():
#          hole_end = ScrewHoles(20e3, origin =primitive.start)
#          if isinstance(primitive, CPW2CPWArc):
#              print('start_angle', primitive.start_angle)
#          hole_end.place(cell, layers[9])    
    #for k in vias:
      #k.place(cell, layers[0])
    
#    print(isinstance(feedline,CPW))
    
#    hole_end = ScrewHoles(20e3, origin =DPoint(25000,-90000))
#    hole_end.place(cell, layers[9])
    
#    hole_arc_1 = ScrewHoles(5e3, origin =DPoint(1000646.987,7612.04674887))
#    hole_arc_1.place(cell, layers[9])
    
#    hole_arc_2 = ScrewHoles(5e3, origin=DPoint(2421166.57537,1406601.51562))
#    hole_arc_2.place(cell, layers[9])
   
   