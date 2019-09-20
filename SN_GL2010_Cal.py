class SN_curve(self):

    def __init__(self,t = 130, *args, **kwargs):
        super().__init__(*args, **kwargs) 
        """
        Basic paramater for SN curve calculation
        """
        # Strength of tension. /MPa
        self.Rm = 360
        self.Rb = Rm*1.06
        # strength of yield. /MPa
        self.Rp = 220
        self.Rs = Rp*1.06
        # Type of the material 1 for spherpodal graphire cast iron; 2 for cast steel.
        self.mat = 1
        # Thickness of the material. /mm
        self.t = t
        # Thickness correction.  ///The thickness is not considered here
        self.sign_tc = 0
        # Surface roughness. 
        self.Rz = 125

    """
    Main function for SN-curve paramater calculation.
    """
    def Fok(self):
        from math import log10, sqrt
        '''
        Influence factor of surface roughness & notch.
        '''
        Rz = self.Rz
        Rb = self.Rb
        # Surface roughness factor.
        #///The details about this factor can be found in 'Leitfaden für eine Betriebsfestigkeitsrechnung'-(Function 8.22).
        Fo = 1-0.22*(log10(Rz)**0.64)*log10(Rb)+0.45*(log10(Rz)**0.53)
        # Notch factor. 
        # ///Which is considered in the FE calculation and Fatigue calculaiton(FEMFAT).
        # ///The details about this factor can be found in 'Leitfaden für eine Betriebsfestigkeitsrechnung'(Chapter 5).
        alphak = 1
        n = 1
        beltk = alphak*n
        # Total influence factor fok.
        # /// The details about this factor can be found in 'Leitfaden für eine Betriebsfestigkeitsrechnung'(Table 8.1).           
        return sqrt(beltk**2-1+1/(Fo**2))

    def Rwk(self):
        """
        Fatigue shtrngth of component
        """
        mat = self.mat
        Rb = self.Rb
        Fok = self.Fok()
        # fatigue strength of specimen
        if mat == 1:
            Rw = 0.27*Rb +100
            M = 0.00035*Rb +0.08
        else:
            Rw = 0.27*Rb +85
            M = 0.00035*Rb +0.05    
        return Rw/Fok