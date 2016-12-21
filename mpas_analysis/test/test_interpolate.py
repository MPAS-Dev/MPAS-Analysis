"""
Tests for the interpolation module

author: Luke Van Roekel
date: 10-24-2016
"""

import numpy as np
from mpas_analysis.shared.interpolation.interpolate import interp_fields, init_tree, lon_lat_to_cartesian
from mpas_analysis.test import TestCase
   
class TestInterp(TestCase):
    def test_lat_to_cartesian(self):
        """
        Test that input lat lon arrays are converted appropriately
    
        Author: Luke Van Roekel
        date: 10-25-2016
        """        
    
        lat_input = np.deg2rad(np.array([-90, 90, 0, 0, 0, 0, 45, 45, 45, 45]))
        lon_input = np.deg2rad(np.array([0, 0, 0, 90, 180, 270, 0, 90, 180, 270]))
        
        x_input = np.cos(lat_input) * np.cos(lon_input)
        y_input = np.cos(lat_input) * np.sin(lon_input)
        z_input = np.sin(lat_input)
        
        x, y, z = lon_lat_to_cartesian(lon_input, lat_input, R=1.0)
        
        self.assertArrayEqual(x, x_input)
        self.assertArrayEqual(y, y_input)
        self.assertArrayEqual(z, z_input)
    
    def test_target_grid_size(self): 
        """
        Test that target lat lon grid is produced with correct size
    
        Author: Luke Van Roekel
        date: 10-25-2016
        """ 
        
        lat_input = np.deg2rad(np.array([0, 0, 0, 0, 45, 45, 45, 45]))
        lon_input = np.deg2rad(np.array([0, 90, 180, 270, 0, 90, 180, 270]))
        
        lon_min = -90
        lon_max = 90
        lat_min = -45
        lat_max = 45
        dLon = 30
        dLat = 30
        
        d, inds, lonTarg, latTarg = init_tree(lon_input, lat_input, lon_min, lon_max, 
                                              lat_min, lat_max, dLon, dLat)
                                              
        nLonExpected = (lon_max - lon_min) / dLon
        nLatExpected = (lat_max - lat_min) / dLat
        
        # LonTarg and LatTarg should be 2-D arrays
        self.assertEqual(len(lonTarg.shape), 2)
        self.assertEqual(len(latTarg.shape), 2)
        
        # LonTarg and latTarg should have expected sizes
        self.assertEqual(lonTarg.shape[0], nLonExpected)
        self.assertEqual(latTarg.shape[0], nLonExpected)
        self.assertEqual(lonTarg.shape[1], nLatExpected)
        self.assertEqual(latTarg.shape[1], nLatExpected)

    def test_target_grid_bounds(self):
        """
        Test that target lat lon grid is produced with correct bounds
    
        Author: Luke Van Roekel
        date: 10-25-2016
        """ 
        
        lat_input = np.deg2rad(np.array([0, 0, 0, 0, 45, 45, 45, 45]))
        lon_input = np.deg2rad(np.array([0, 90, 180, 270, 0, 90, 180, 270]))
        
        lon_min = -180
        lon_max = 180
        lat_min = -50
        lat_max = 50
        dLon = 20
        dLat = 20
        
        d, inds, lonTarg, latTarg = init_tree(lon_input, lat_input, lon_min, lon_max, 
                                              lat_min, lat_max, dLon, dLat)  
                             
        # LonTarg should respect bounds defined by user
        self.assertLessThan(lonTarg.max(), lon_max)
        self.assertGreaterThan(lonTarg.min(), lon_min)
        
        self.assertLessThan(latTarg.max(), lat_max)
        self.assertGreaterThan(latTarg.min(), lat_min)
        
    def test_interp(self):
        """
        Test that nearest neighbor interpolation works as expected
    
        Author: Luke Van Roekel
        date: 10-25-2016
        """         
        lat_input = np.deg2rad(np.array([-45, -45, 45, 45]))
        lon_input = np.deg2rad(np.array([-90, 90, -90, 90]))
        vals_input = np.array([1, 2, 3, 4])
        
        lon_min = -90
        lon_max = 90
        lat_min = -45
        lat_max = 45
        dLon = 30
        dLat = 30
        
        d, inds, lonTarg, latTarg = init_tree(lon_input, lat_input, lon_min, lon_max, 
                                              lat_min, lat_max, dLon, dLat)  
                                              
        vals_output = interp_fields(vals_input, d, inds, lonTarg)
        
        #Test a few ramdom spots for nearest neighbor
        #One exact spot
        self.assertEqual(vals_output[0, 0], vals_input[0])
        
        # Close point to val_input(0)
        self.assertEqual(vals_output[2, 0], vals_input[0])
        
        # Close to point val_input(1)
        self.assertEqual(vals_output[4, 0], vals_input[1])
        
        #Close to val_input(2)
        self.assertEqual(vals_output[2, 2], vals_input[2])
        
        #Close to val_input(3)
        self.assertEqual(vals_output[4, 2], vals_input[3])
    