from math import sqrt, cos, sin, atan2, pi, copysign, tan

import pya
import numpy as np
from pya import Point, DPoint, DVector, DSimplePolygon, SimplePolygon, \
    DPolygon, Polygon, Region
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from typing import Union, List
from collections import OrderedDict
import itertools
import copy

from classLib.baseClasses import ElementBase, ComplexBase
from classLib.bridgedCoplanars import BridgedCPW, BridgedCPWArc


class CirclePolygon():
  
    def __init__(self, radius, number_of_points=10):
        self.radius = radius
        self.number_of_points = number_of_points
    
    def get_circle(self):
        radius = self.radius
        nr_points = self.number_of_points
        angles = np.linspace(0,2*np.pi,nr_points+1)[0:-1]
        points = []

        for ind,angle in enumerate(angles):
            points.append(pya.Point(radius*np.cos(angle),radius*np.sin(angle)))
            circle = pya.SimplePolygon(points)
      
        return circle


class ScrewHoles(ElementBase):

    def __init__(self, radius=None, origin=DPoint(0, 0), trans_in=None, region_id="default"):
       
        self.radius = radius
        super().__init__(
            origin=origin,
            trans_in=trans_in,
            region_id=region_id
        )

    def init_regions(self):
        circle = CirclePolygon(self.radius).get_circle()
        self.empty_region.insert(circle)
  
       # self.empty_region.size(1, 0, 0)


class MetalRounds(ElementBase):

    def __init__(self, radius=None, origin=DPoint(0, 0), trans_in=None, region_id="default"):
        self.radius = radius
        super().__init__(
            origin=origin,
            trans_in=trans_in,
            region_id=region_id
        )

    def init_regions(self):
        circle = CirclePolygon(self.radius).get_circle()
        self.metal_region.insert(circle)

    # self.empty_region.size(1, 0, 0)
       
       
class SampleHoles(ElementBase):

    def __init__(self, dx = None, dy = None, step_forward = None, origin = DPoint(0,0), gap_from_step = None,
               trans_in = None, ports_num = 40, ports = [], cpw_params = None, rounds_radius = None, via_r=None,
                 region_id = "default"):
        """
        CPW_PARAMS['g'] - gap width (two gaps ofc)
        CPW_PARAMS['s'] - central line width
        CPW_PARAMS['w'] - shortest distance from the gap edge to the via edge
        step_forward - step for metallization from chip hole
        gap_from_step - gap from metallization edge (for connection pads)

        via_width (distance between centers of circles) from via class should be parametrized as:
        via_width = via_r+cpw_params['w']+cpw_params['g']+cpw_params['s']+cpw_params['g']+cpw_params['w']+via_r =
        = 2*via_r + 2*cpw_params['w'] + 2*cpw_params['g'] + cpw_params['s']

        ports_num is hardcoded for 40 pin holder

        ports list contains DPoints for the  CPWOPENED start point
        """
        self.dx = dx
        self.dy = dy
        self.gap_from_step = gap_from_step
        self.step_forward = step_forward
        self.via_r = via_r
        self.ports_num = ports_num
        self.ports = ports
        self.rounds_radius = rounds_radius
        self.cpw_params = cpw_params
        self.empty_box_coords = []
        self.empty_box_coords_steped_forward = []

        super().__init__(origin = origin, trans_in = trans_in, region_id = region_id)

    def init_regions(self):

        _empty_box_coords = [DPoint(-self.dx / 2, self.dy / 2),
                               DPoint(self.dx / 2, self.dy / 2),
                               DPoint(self.dx / 2, -self.dy / 2),
                               DPoint(-self.dx / 2, - self.dy / 2)]

        _empty_box_coords_steped_forward = [DPoint(-self.dx / 2 - self.step_forward, self.dy / 2 + self.step_forward),
                               DPoint(self.dx / 2 + self.step_forward, self.dy / 2 + self.step_forward),
                               DPoint(self.dx / 2 + self.step_forward, -self.dy / 2 - self.step_forward),
                               DPoint(-self.dx / 2 - self.step_forward, - self.dy / 2 - self.step_forward)]

        self.empty_box_coords_steped_forward = _empty_box_coords_steped_forward
        self.empty_box_coords = _empty_box_coords

        _empty_box = DPolygon(_empty_box_coords)

        _empty_box_stepped_forward = DPolygon(_empty_box_coords_steped_forward)

        _geometries = [_empty_box]+[_empty_box_stepped_forward]

        for geometry in _geometries:
            self.empty_region.insert(geometry)

        ports_num = self.ports_num
        ports = self.ports
        cpw_params = self.cpw_params
        via_r = self.via_r
        self.step_forward = self.step_forward + self.gap_from_step
        offset = (self.dx + 2 * self.step_forward) / 11  #cpw_params['s']*7 + cpw_params['w']*2+cpw_params['g']*2 + via_r*4

        for port in range(1, ports_num+1):
            if port<= 10:
                ports.append(
                    DPoint(-self.dx / 2 - self.step_forward + port * offset, self.dy / 2 + self.step_forward))
            elif 11<=port<=20:
                ports.append(
                    DPoint(self.dx / 2 + self.step_forward, self.dy / 2 + self.step_forward-(port-10)*offset))
            elif 21<=port<=30:
                ports.append(
                    DPoint(self.dx / 2 + self.step_forward - (port-20) * offset, -self.dy / 2 - self.step_forward))
            elif 31<=port<=40:
                ports.append(
                    DPoint(-self.dx / 2 - self.step_forward, - self.dy / 2 - self.step_forward + (port - 30) * offset)
                )

            self.ports = ports

