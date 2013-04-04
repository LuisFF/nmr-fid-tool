package uk.ac.ebi.nmr.fid.tools;

import uk.ac.ebi.nmr.fid.Proc;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 18/02/2013
 * Time: 05:29
 * To change this template use File | Settings | File Templates.
 */
public class PhaseCorrectionTool {

    double [] spectrum;
    Proc processing;
    double [] data;

    public PhaseCorrectionTool(double[] spectrum, Proc processing) {
        this.spectrum=spectrum;
        this.processing=processing;
        this.data = new double[processing.getTdEffective()/2+1];
    }

    public double [] zeroOrderPhasing (double angle){
        for(int i =0; i < processing.getTdEffective()-1 ; i+=2){
            data[i/2]=spectrum[i]*Math.cos(angle)-spectrum[i+1]*Math.sin(angle);
        }
        return data;
    }

    public double [] firstOrderPhasing (double angle){
        for(int i=0; i< processing.getTdEffective()-1; i+=2){
            data[i/2]=spectrum[i]*Math.cos(i/2*angle/(processing.getTdEffective()/2))
                    +spectrum[i+1]*Math.sin(i/2*angle/(processing.getTdEffective()/2));
        }
        return data;
    }




}
