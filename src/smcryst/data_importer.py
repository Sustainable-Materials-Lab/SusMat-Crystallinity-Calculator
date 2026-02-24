#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data import handler for WAXS crystallinity calculator

This module contains the DataImporter class that handles the complex data import
and parsing logic for different phase configurations and file formats.
"""
import copy
import numpy as np


class DataImporter:
    """
    Handles data import and parsing for XRD/WAXS data with various phase configurations.
    
    This class encapsulates the complex logic for parsing data files with different
    phase combinations (Cellulose I/II, iPP, PCL, etc.) and various file formats.
    """
    
    def __init__(self, args):
        """
        Initialize the DataImporter with command line arguments.
        
        Args:
            args: Parsed command line arguments containing configuration options
        """
        self.args = args
        self.data = None
        self.raw_text = None
        
        # Initialize phase flags based on file content or args
        self.sqrt = False
        self.iPP = False
        self.giPP = False
        self.PCL = False
        self.cel2 = args.cel2 if hasattr(args, 'cel2') else False
        self.jeffamine = False
        self.cel2amorph = False
        self.background = False
        self.expback = False
        self.xye = False
        
    def load_data(self, input_data):
        """
        Load data from file or raw text and detect phase configurations.
        
        Args:
            input_data: Either a path to the input data file (str) or raw text lines (list)
            
        Returns:
            dict: Dictionary containing all parsed data arrays
        """
        if isinstance(input_data, str):
            # Load data from file path
            self.data = np.genfromtxt(input_data, delimiter=',', skip_header=2,
                                    invalid_raise=False, unpack=False)
            # Read the header to detect phases
            self._detect_phases(input_data)
        else:
            # Load data from raw text lines
            self.data = np.genfromtxt(input_data, delimiter=',', skip_header=2,
                                    invalid_raise=False, unpack=False)
            # Set raw_text for phase detection
            self.raw_text = input_data
            # Detect phases based on raw text
            self._detect_phases_from_text()
        
        # Parse the data based on detected phases and configurations
        parsed_data = self._parse_data()
        
        # Apply exposure correction
        self._apply_exposure_correction(parsed_data)
        
        # Apply clipping if requested
        if self.args.clip:
            parsed_data = self._apply_clipping(parsed_data)
            
        # Apply linear subtraction if requested
        if self.args.linsub:
            parsed_data = self._apply_linear_subtraction(parsed_data)
            
        # Apply background handling
        if self.args.keepbkg:
            parsed_data = self._apply_background_restoration(parsed_data)
            
        return parsed_data
    
    def _detect_phases(self, filename):
        """
        Detect which phases are present based on file header content.
        
        Args:
            filename: Path to the input data file
        """
        with open(filename, encoding='utf-8') as f:
            self.raw_text = f.readlines()

        if "Sqrt(y)" in self.raw_text[0]:
            self.sqrt = True
            data2 = np.square(self.data)
            data2 = np.delete(data2, 0, 1)
            self.data = np.insert(data2, 0, self.data[:, 0], 1)

        # Detect phases from header
        header_line = self.raw_text[1] if len(self.raw_text) > 1 else ""
        
        self.iPP = "Alpha i-PP" in header_line
        self.giPP = "Gamma i-PP" in header_line
        self.PCL = "PCL" in header_line
        self.jeffamine = "Jeffamine ED2003" in header_line
        self.cel2amorph = "Amorphous,Background" in header_line
        self.background = "Bkg" in header_line
        self.expback = "Background,Amorphous" in header_line
        
        # Detect if file has errors included (xye format)
        self.xye = "SigmaYobs" in header_line
        
        # Cellulose II detection (can be overridden by args)
        if "Cellulose II" in header_line:
            self.cel2 = True
        elif hasattr(self.args, 'cel2') and self.args.cel2:
            self.cel2 = True
    
    def _detect_phases_from_text(self):
        """
        Detect which phases are present based on raw text lines.
        Used when input_data is a list of text lines instead of a filename.
        """
        if "Sqrt(y)" in self.raw_text[0]:
            self.sqrt = True
            data2 = np.square(self.data)
            data2 = np.delete(data2, 0, 1)
            self.data = np.insert(data2, 0, self.data[:, 0], 1)

        # Detect phases from header
        header_line = self.raw_text[1] if len(self.raw_text) > 1 else ""
        
        self.iPP = "Alpha i-PP" in header_line
        self.giPP = "Gamma i-PP" in header_line
        self.PCL = "PCL" in header_line
        self.jeffamine = "Jeffamine ED2003" in header_line
        self.cel2amorph = "Amorphous,Background" in header_line
        self.background = "Bkg" in header_line
        self.expback = "Background,Amorphous" in header_line
        
        # Detect if file has errors included (xye format)
        self.xye = "SigmaYobs" in header_line
        
        # Cellulose II detection (can be overridden by args)
        if "Cellulose II" in header_line:
            self.cel2 = True
        elif hasattr(self.args, 'cel2') and self.args.cel2:
            self.cel2 = True
    
    def _parse_data(self):
        """
        Parse the data array based on detected phases and configurations.
        
        Returns:
            dict: Dictionary containing all parsed data arrays
        """
        data = self.data.T
        result = {}
        
        if self.args.peaks:
            result = self._parse_peaks_data(data)
        elif self.background:
            result = self._parse_background_data(data)
        elif self.expback:
            result = self._parse_expback_data(data)
        else:
            result = self._parse_standard_data(data)
            
        return result
    
    def _parse_peaks_data(self, data):
        """Parse data when amorphous phase is modelled as xo_Is peaks."""
        result = {}
        
        if self.xye:
            if self.cel2:
                result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['cel2'] = data
                result['amorph'] = result['prf'] - result['cel1'] - result['cel2']
            elif self.iPP:
                if self.giPP:
                    result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['gipp'] = data
                    result['amorph'] = result['prf'] - result['cel1'] - result['gipp'] - result['iPP']
                else:
                    result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['iPP'] = data
                    result['amorph'] = result['prf'] - result['cel1'] - result['iPP']
            elif self.PCL:
                result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['PCL'] = data
                result['amorph'] = result['prf'] - result['cel1'] - result['PCL']
            else:
                result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'] = data
                result['amorph'] = result['prf'] - result['cel1']
        else:
            if self.cel2:
                result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['cel2'] = data
                result['amorph'] = result['prf'] - result['cel1'] - result['cel2']
            elif self.iPP:
                if self.giPP:
                    result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['gipp'] = data
                    result['amorph'] = result['prf'] - result['cel1'] - result['gipp'] - result['iPP']
                else:
                    result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['iPP'] = data
                    result['amorph'] = result['prf'] - result['cel1'] - result['iPP']
            elif self.PCL:
                result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['PCL'] = data
                result['amorph'] = result['prf'] - result['cel1'] - result['PCL']
            else:
                result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'] = data
                result['amorph'] = result['prf'] - result['cel1']
                
        # Initialize bck as zeros
        result['bck'] = np.zeros_like(result['ttheta'])
        return result
    
    def _parse_background_data(self, data):
        """Parse data when background is defined."""
        result = {}
        
        if hasattr(self.args, 'bckamorph') and self.args.bckamorph:
            # Background models amorphous content
            if self.xye:
                if self.cel2:
                    result['ttheta'], result['raw'], result['_'], result['bck'], result['prf'], result['diff'], result['cel1'], result['cel2'] = data
                    result['cel2'] -= result['bck']
                elif self.iPP:
                    if self.giPP:
                        result['ttheta'], result['raw'], result['_'], result['bck'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['gipp'] = data
                        result['iPP'] -= result['bck']
                        result['gipp'] -= result['bck']
                    else:
                        result['ttheta'], result['raw'], result['_'], result['bck'], result['prf'], result['diff'], result['cel1'], result['iPP'] = data
                        result['iPP'] -= result['bck']
                elif self.PCL:
                    result['ttheta'], result['raw'], result['_'], result['bck'], result['prf'], result['diff'], result['cel1'], result['PCL'] = data
                    result['PCL'] -= result['bck']
                else:
                    result['ttheta'], result['raw'], result['_'], result['bck'], result['prf'], result['diff'], result['cel1'] = data
            else:
                if self.cel2:
                    result['ttheta'], result['raw'], result['bck'], result['prf'], result['diff'], result['cel1'], result['cel2'] = data
                    result['cel2'] -= result['bck']
                elif self.iPP:
                    if self.giPP:
                        result['ttheta'], result['raw'], result['bck'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['gipp'] = data
                        result['iPP'] -= result['bck']
                        result['gipp'] -= result['bck']
                    else:
                        result['ttheta'], result['raw'], result['bck'], result['prf'], result['diff'], result['cel1'], result['iPP'] = data
                        result['iPP'] -= result['bck']
                elif self.PCL:
                    result['ttheta'], result['raw'], result['bck'], result['prf'], result['diff'], result['cel1'], result['PCL'] = data
                    result['PCL'] -= result['bck']
                else:
                    result['ttheta'], result['raw'], result['bck'], result['prf'], result['diff'], result['cel1'] = data

            result['cel1'] -= result['bck']
            result['amorph'] = copy.copy(result['bck'])
            result['bck'] *= 0
        else:
            # Standard background subtraction
            if self.xye:
                if self.cel2:
                    result['ttheta'], result['raw'], result['_'], result['bck'], result['prf'], result['diff'], result['cel1'], result['cel2'], result['amorph'] = data
                elif self.iPP:
                    if self.giPP:
                        result['ttheta'], result['raw'], result['_'], result['bck'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['gipp'], result['amorph'] = data
                    else:
                        result['ttheta'], result['raw'], result['_'], result['bck'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['amorph'] = data
                elif self.PCL:
                    result['ttheta'], result['raw'], result['_'], result['bck'], result['prf'], result['diff'], result['cel1'], result['PCL'], result['amorph'] = data
                elif self.jeffamine:
                    result['ttheta'], result['raw'], result['_'], result['bck'], result['prf'], result['diff'], result['cel1'], result['jeffamine'], result['amorph'] = data
                else:
                    result['ttheta'], result['raw'], result['_'], result['bck'], result['prf'], result['diff'], result['cel1'], result['amorph'] = data
            else:
                if self.cel2:
                    result['ttheta'], result['raw'], result['bck'], result['prf'], result['diff'], result['cel1'], result['cel2'], result['amorph'] = data
                elif self.iPP:
                    if self.giPP:
                        result['ttheta'], result['raw'], result['bck'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['gipp'], result['amorph'] = data
                    else:
                        result['ttheta'], result['raw'], result['bck'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['amorph'] = data
                elif self.PCL:
                    result['ttheta'], result['raw'], result['bck'], result['prf'], result['diff'], result['cel1'], result['PCL'], result['amorph'] = data
                else:
                    result['ttheta'], result['raw'], result['bck'], result['prf'], result['diff'], result['cel1'], result['amorph'] = data
            
            # Apply background subtraction
            self._subtract_background(result)
            
        return result
    
    def _parse_expback_data(self, data):
        """Parse data when experimental background mode is enabled."""
        result = {}
        
        if self.xye:
            if self.cel2:
                result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['cel2'], result['bck'], result['amorph'] = data
            elif self.iPP:
                if self.giPP:
                    result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['gipp'], result['bck'], result['amorph'] = data
                else:
                    result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['bck'], result['amorph'] = data
            elif self.PCL:
                result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['PCL'], result['bck'], result['amorph'] = data
            else:
                result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['bck'], result['amorph'] = data
        else:
            if self.cel2:
                result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['cel2'], result['bck'], result['amorph'] = data
            elif self.iPP:
                if self.giPP:
                    result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['gipp'], result['bck'], result['amorph'] = data
                else:
                    result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['bck'], result['amorph'] = data
            elif self.PCL:
                result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['PCL'], result['bck'], result['amorph'] = data
            else:
                result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['bck'], result['amorph'] = data
        
        # Apply background subtraction
        self._subtract_background(result)
        
        return result
    
    def _parse_standard_data(self, data):
        """Parse data in standard mode."""
        result = {}
        
        if self.xye:
            if self.cel2:
                if self.cel2amorph:
                    result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['cel2'], result['amorph'] = data
                else:
                    result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['cel2'], result['amorph'] = data
            elif self.iPP:
                if self.giPP:
                    if self.cel2amorph:
                        result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['amorph'], result['gipp'] = data
                    else:
                        result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['amorph'], result['gipp'] = data
                else:
                    if self.cel2amorph:
                        result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['amorph'] = data
                    else:
                        result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['amorph'] = data
            elif self.PCL:
                if self.cel2amorph:
                    result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['PCL'], result['amorph'] = data
                else:
                    result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['PCL'], result['amorph'] = data
            else:
                if self.cel2amorph:
                    result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['amorph'] = data
                else:
                    result['ttheta'], result['raw'], result['_'], result['prf'], result['diff'], result['cel1'], result['amorph'] = data
        else:
            if self.cel2:
                if self.cel2amorph:
                    result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['cel2'], result['amorph'] = data
                else:
                    result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['cel2'], result['amorph'] = data
            elif self.iPP:
                if self.giPP:
                    if self.cel2amorph:
                        result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['amorph'], result['gipp'] = data
                    else:
                        result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['amorph'], result['gipp'] = data
                else:
                    if self.cel2amorph:
                        result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['amorph'] = data
                    else:
                        result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['iPP'], result['amorph'] = data
            elif self.PCL:
                if self.cel2amorph:
                    result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['PCL'], result['amorph'] = data
                else:
                    result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['PCL'], result['amorph'] = data
            else:
                if self.cel2amorph:
                    result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['amorph'] = data
                else:
                    result['ttheta'], result['raw'], result['prf'], result['diff'], result['cel1'], result['amorph'] = data
        
        # Initialize bck as zeros if not present
        if 'bck' not in result:
            result['bck'] = np.zeros_like(result['ttheta'])
            
        return result
    
    def _subtract_background(self, result):
        """Apply background subtraction to all relevant arrays."""
        bck = result['bck']
        
        result['prf'] -= bck
        result['raw'] -= bck
        result['amorph'] -= bck
        result['cel1'] -= bck
        
        if 'cel2' in result:
            result['cel2'] -= bck
        if 'iPP' in result:
            result['iPP'] -= bck  
        if 'gipp' in result:
            result['gipp'] -= bck
        if 'PCL' in result:
            result['PCL'] -= bck
        if 'jeffamine' in result:
            result['jeffamine'] -= bck
    
    def _apply_exposure_correction(self, result):
        """Apply exposure time correction to all intensity arrays."""
        exposure = self.args.exposure
        
        result['raw'] *= exposure
        result['prf'] *= exposure
        result['diff'] *= exposure
        result['cel1'] *= exposure
        result['amorph'] *= exposure
        
        if 'cel2' in result:
            result['cel2'] *= exposure
        if 'iPP' in result:
            result['iPP'] *= exposure
        if 'gipp' in result:
            result['gipp'] *= exposure
        if 'PCL' in result:
            result['PCL'] *= exposure
        if 'jeffamine' in result:
            result['jeffamine'] *= exposure
    
    def _apply_clipping(self, result):
        """Apply 2theta range clipping."""
        ttheta = result['ttheta']
        startval = np.where(ttheta > self.args.xmin)[0][0]
        endval = np.where(ttheta > self.args.xmax)[0][0]
        
        # Apply clipping to all arrays
        for key, array in result.items():
            if isinstance(array, np.ndarray) and len(array) == len(ttheta):
                result[key] = array[startval:endval]
                
        return result
    
    def _apply_linear_subtraction(self, result):
        """Remove linear portion of amorphous scattering and add to background."""
        # Remove minimum values from each component
        minamorph = min(result['amorph'])
        result['amorph'] -= minamorph
        
        mincel1 = min(result['cel1'])
        result['cel1'] -= mincel1
        
        adjustments = [minamorph, mincel1]
        
        if 'cel2' in result:
            mincel2 = min(result['cel2'])
            result['cel2'] -= mincel2
            adjustments.append(mincel2)
            
        if 'iPP' in result:
            miniPP = min(result['iPP'])
            result['iPP'] -= miniPP
            adjustments.append(miniPP)
            
        if 'gipp' in result:
            mingipp = min(result['gipp'])
            result['gipp'] -= mingipp
            adjustments.append(mingipp)
            
        if 'PCL' in result:
            minpcl = min(result['PCL'])
            result['PCL'] -= minpcl
            adjustments.append(minpcl)
            
        if 'jeffamine' in result:
            minjeffamine = min(result['jeffamine'])
            result['jeffamine'] -= minjeffamine
            adjustments.append(minjeffamine)
        
        # Apply total adjustment to profile and raw data
        total_adjustment = sum(adjustments)
        result['prf'] -= total_adjustment
        result['raw'] -= total_adjustment
        
        return result
    
    def _apply_background_restoration(self, result):
        """Add background back to all components for plotting."""
        bck = result.get('bck', np.zeros_like(result['ttheta']))
        
        result['raw'] += bck
        result['cel1'] += bck
        result['prf'] += bck
        result['amorph'] += bck
        
        if 'cel2' in result:
            result['cel2'] += bck
        if 'iPP' in result:
            result['iPP'] += bck
        if 'gipp' in result:
            result['gipp'] += bck
        if 'PCL' in result:
            result['PCL'] += bck
        if 'jeffamine' in result:
            result['jeffamine'] += bck
            
        return result