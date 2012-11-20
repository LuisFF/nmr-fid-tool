/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.ebi.nmr;

import java.util.Collections;
import java.util.List;

/**
 *
 * @author ldpf
 */
public class SpectrumGenerator {

    private static List<Double> shifts;
    private static List<Double> intensities;
    private static Double[] x_curve;
    private static Double[] y_curve;
    private static Double[] x_peaks;
    private static Double[] y_peaks;
    private static Double scaleX;
    private static Double scaleY;
    private static Double shiftX;
    private static Double min_shift;
    private static Double max_shift;

    public SpectrumGenerator(List<Double> shifts, List<Double> intensities) {
        this.shifts = shifts;
        this.intensities = intensities;
        this.scaleX=1.0;
        this.scaleY=1.0;
        //shifts the shifts x units to the left
        this.shiftX=0.0;
        this.min_shift = Collections.min(shifts);
        this.max_shift = Collections.max(shifts);
        
    }

    public void computeCurvePoints() throws Exception {
        int i;
        int j;
	double lineshape= 0.0002;
        double line_width = 0.2;
        int point_size = 0;
//	GdkColor line_color = winData->line_color;
//	GdkColor point_color = winData->point_color;
        int n = 0;
        
        double xx = min_shift;
        double halfWidth = Math.abs(max_shift - min_shift) / 30;
        double h0 = halfWidth / 20;

//	if(dataCurve->x && dataCurve->y)
//	{
//		line_width = dataCurve->line_width;
//		point_size = dataCurve->point_size;
//		line_color = dataCurve->line_color;
//		point_color = dataCurve->point_color;
//	}
//
//	if(winData->convType==GABEDIT_CONV_TYPE_LORENTZ)
//	       lineshape = lorentzianLineshape;
//	else
//	       lineshape = gaussianLineshape;
//
//
//	
//	if(dataCurve->size>0 && winData->size)
//	{
//		if(dataCurve->x) g_free(dataCurve->x);
//		if(dataCurve->y) g_free(dataCurve->y);
//		dataCurve->x = NULL;
//		dataCurve->y = NULL;
//	}

//	dataCurve->size=0;
        int dataCurve_size = 0;
//        xx = winData->xmin;
        // step size at which we will progress along the x axis
        h0 = halfWidth / 10;


        if (shifts.isEmpty() || shifts.size() < 1 || intensities.isEmpty() || intensities.size() < 1) {
            throw new Exception(SpectrumGenerator.class + ": No shifts and intensities defined");
        }

        n = (int) (((max_shift - min_shift) / h0 + 0.5) + shifts.size());
        
        Double[] x_tmp =null;
        if (n > 0) {
            x_tmp = new Double[n];
        }
        // define the shift center?
        if (shifts.size() > 0 && n > 0) {
            do {
                double dmin = 0.0;
                double d = 0.0;
                int jmin = 0;
                for (j = 0; j < shifts.size(); j++) {
                    double center = shifts.get(j) *  scaleX + shiftX;

                    d = Math.abs(xx - center);
                    if (d < dmin || j == 0) {
                        jmin = j;
                        dmin = d;
                    }
                }
                x_tmp[dataCurve_size] = xx;
                if (dmin < h0) {
                    if (xx < shifts.get(jmin) * scaleX + shiftX) {
                        xx = shifts.get(jmin) * scaleX + shiftX;
                        x_tmp[dataCurve_size] = xx;
                        xx += h0 + 1e-8;
                    } else {
                        xx =  shifts.get(jmin) * scaleX + shiftX;
                        x_tmp[dataCurve_size] = xx;
                        xx += h0 + 1e-8;
                    }
                } else {
                    if (dmin > 5 * halfWidth) {
                        xx += h0 + dmin / 5;
                    } else {
                        xx += h0;
                    }
                }

                dataCurve_size++;
            } while (xx < max_shift && dataCurve_size < n);
        }
        
        // rescale x and define size of y
        if (dataCurve_size > 0) {
            x_curve = new Double[dataCurve_size];
            y_curve = new Double[dataCurve_size];
        }
//        System.out.println(dataCurve_size+" "+x_tmp.length);
        // pass all the points that we are interested on to the final x
        // not that the size of the x_tmp should be bigger than x
        for (int k = 0; k < dataCurve_size; k++) {
            x_curve[k]=x_tmp[k];
        }
        
        if (shifts.size() > 0) {
            for (i = 0; i < dataCurve_size; i++) {
                Double yy = 0.0;
                for (j = 0; j < shifts.size(); j++) {
                    double center =  shifts.get(j) * scaleX + shiftX;
                    double rel_offset = (x_curve[i] - center) / halfWidth;
                    yy += intensities.get(j) * lineshape * scaleY;
                }
                y_curve[i] = yy;
            }
        }

//        sprintf(dataCurve - > point_str, "+");
//        dataCurve - > point_size = point_size;
//        dataCurve - > line_width = line_width;
//        dataCurve - > point_color = point_color;
//
//        dataCurve - > line_color = line_color;
//
//        dataCurve - > line_style = winData - > line_style;
    }
    
    public static void computePeaks() throws Exception{
        int loop;
         
//	gint line_width = winData->line_width;
//	gint point_size = winData->point_size;
//	GdkColor line_color = winData->line_color;
//	GdkColor point_color = winData->point_color;
	
//	line_color.red *=0.9; 
//	line_color.green *=0.9; 
//	line_color.blue *=0.9; 
//	point_color.red *=0.9; 
//	point_color.green *=0.9; 
//	point_color.blue *=0.9; 
//	if(dataPeaks->x && dataPeaks->y)
//	{
//		line_width = dataPeaks->line_width;
//		point_size = dataPeaks->point_size;
//		line_color = dataPeaks->line_color;
//		point_color = dataPeaks->point_color;
//	}
        

	int dataPeaks_size=3*shifts.size()+2;
        if (dataPeaks_size<1) {
            throw new Exception(SpectrumGenerator.class + ": No shifts and intensities defined");
        }
	
        x_peaks = new Double[dataPeaks_size];
        y_peaks = new Double[dataPeaks_size];
		
	
	
     
	x_peaks[0]=min_shift;
	/* dataPeaks->y[0]=winData->ymin;*/
	y_peaks[0]=0.0;
	x_peaks[dataPeaks_size-1]=max_shift;
	/* dataPeaks->y[dataPeaks->size-1]=winData->ymin;*/
	y_peaks[dataPeaks_size-1]=0.0;
	for (loop=0; loop<shifts.size(); loop++){
		int iold = loop*3+1;
		double xx = shifts.get(loop)*scaleX+shiftX;
		x_peaks[iold]=xx;
		/* dataPeaks->y[iold]=winData->ymin;*/
		y_peaks[iold]=0.0;

		x_peaks[iold+1]=xx;
		y_peaks[iold+1]=intensities.get(loop)*scaleY;

		x_peaks[iold+2]=xx;
		/* dataPeaks->y[iold+2]=winData->ymin;*/
		y_peaks[iold+2]=0.0;
	}
//
//	sprintf(dataPeaks->point_str,"+");
//	dataPeaks->point_size=point_size;
//	dataPeaks->line_width=line_width;
//	dataPeaks->point_color=point_color; 
//
//	dataPeaks->line_color=line_color; 
//	dataPeaks->line_style=winData->line_style; 
    }
    public static void setYmax2One() throws Exception{
        double ymax = 0;
	int loop;
	if(x_curve.length<1) 
            throw new Exception("No curve points available");
        
	ymax = y_curve[0];
        //??? I only get a horizontal line...
	for (loop=1; loop<y_curve.length; loop++)
		if(ymax<y_curve[loop]) ymax = y_curve[loop];
	if(ymax!=0)
	for (loop=0; loop<y_curve.length; loop++)
		y_curve[loop] /= ymax;
	if(ymax!=0)
	for (loop=0; loop<y_peaks.length; loop++)
		y_peaks[loop] /= ymax;
    }

    public static Double[] getX_curve() {
        return x_curve;
    }

    public static Double[] getY_curve() {
        return y_curve;
    }
    public static Double[] getX_peaks() {
        return x_peaks;
    }

    public static Double[] getY_peaks() {
        return y_peaks;
    }
    
}
