package uk.ac.ebi.nmr.fid.tools;

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;

/**
 * The ApodizationTool allows to apply various methods to improve the signal to noise value of an fid.
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 14/01/2013
 * Time: 13:54
 * To change this template use File | Settings | File Templates.
 */
public class ApodizationTool {
    private static double[] spectrum;
    private Acqu.AcquisitionMode acquisitionMode= Acqu.AcquisitionMode.SEQUENTIAL;
    private Proc processing;

    /**
     *
     * @param spectrum
     * @param processing
     */
    public ApodizationTool(double[] spectrum, Proc processing) {
        this.spectrum=spectrum;
        this.processing=processing;
    }
    // TODO perhaps change this to accept an Acqu object, to make aparent that the acquisition mode has to be defined
    public ApodizationTool(double[] spectrum, Acqu.AcquisitionMode acquisitionMode, Proc processing) {
        this.spectrum=spectrum;
        this.processing=processing;
        this.acquisitionMode=acquisitionMode;
    }

    public static double[] getSpectrum() {
        return spectrum;
    }

    /**
     * performs the exponential apodization of the fid using the function F(x,i)= x.exp(-i*dw*lb)
     * where i*dw gives the time coordinate in (s) and the line broadening is in Hz.
     * The cuteNMR implementation is W(x,i) = x.exp(-(dw*i) * lb * pi)
     * @return
     */
    public double [] exponential(){
        double data[] = new double[processing.getTdEffective()];
        if (acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i++)
//                spectrum[i]*= Math.exp(-i*processing.getDwellTime()*Math.PI*processing.getLineBroadning());
                data[i] = spectrum[i]* Math.exp(-i*processing.getDwellTime()*processing.getLineBroadning());
        } else { // simultaneous data
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i+=2){
//                double factor = Math.exp(-i*processing.getDwellTime()*Math.PI*processing.getLineBroadning());
                double factor = Math.exp(-i*processing.getDwellTime()*processing.getLineBroadning());
                data[i]= spectrum[i]* factor;
                data[i+1] =spectrum[i+1]* factor;
//                spectrum[i+1]*= Math.exp(-i*processing.getDwellTime()*Math.PI*processing.getLineBroadning());
            }
        }
        return data;
    }


    public void lorentzGaus(){
        lorentzGaus(0);
    }

    public void lorentzGaus(double lbLorentzToGauss){

        double acquisitionTime = processing.getDwellTime()*processing.getTdEffective();
        double a = Math.PI*lbLorentzToGauss;
        double b = -a/(2.0*processing.getGbFactor()*acquisitionTime);
        double time;

        if (acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i++){
                time=i*processing.getDwellTime();
                spectrum[i]*=Math.exp(-a*time -b*time*time);
            }
        } else {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i+=2){
                time=i*processing.getDwellTime();
                spectrum[i]*= Math.exp(-a*time -b*time*time);
                spectrum[i+1]*= Math.exp(-a*time -b*time*time);
            }
        }

    }

    public void tarf(){
        traf(0);
    }

    public void traf(double lbTraf){
        double acquisitionTime = processing.getDwellTime()*processing.getTdEffective();
        double time, e , eps;

        if (acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i++){
                time=i*processing.getDwellTime();
                e = Math.exp(-time*lbTraf*Math.PI);
                eps = Math.exp((time-acquisitionTime) * lbTraf * Math.PI);
                spectrum[i]*=e/(e*e+eps*eps);
            }
        } else {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i+=2){
                time=i*processing.getDwellTime();
                e = Math.exp(-time*lbTraf*Math.PI);
                eps = Math.exp((time-acquisitionTime) * lbTraf * Math.PI);
                spectrum[i]*=e/(e*e+eps*eps);
                spectrum[i+1]*=e/(e*e+eps*eps);

            }
        }
    }

    public void trafs (){
        trafs(0);
    }

    public void trafs (double lbTraf){
        double acquisitionTime = processing.getDwellTime()*processing.getTdEffective();
        double time, e , eps;

        if (acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i++){
                time=i*processing.getDwellTime();
                e = Math.exp(-time*lbTraf*Math.PI);
                eps = Math.exp((time-acquisitionTime) * lbTraf * Math.PI);
                spectrum[i]*=(e*e*(e+eps))/(e*e*e+eps*eps*eps);
            }
        } else {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i+=2){
                time=i*processing.getDwellTime();
                e = Math.exp(-time*lbTraf*Math.PI);
                eps = Math.exp((time-acquisitionTime) * lbTraf * Math.PI);
                spectrum[i]*=(e*e*(e+eps))/(e*e*e+eps*eps*eps);
                spectrum[i+1]*=(e*e*(e+eps))/(e*e*e+eps*eps*eps);
            }
        }
    }

    public void sine() {
        double factor;
        // ssbSine is defined as 180/ssb
        //
        double offset = (180.0 - processing.getSsbSine())/180.0;
        // check the loop conditions in the original code...
        for (int i =0 ;
                i< (processing.getTdEffective()-1) && i >= ((processing.getTdEffective()-1)-spectrum.length) ;
                i+=2){
            factor=Math.sin(i/processing.getTdEffective()*Math.PI*offset);
            spectrum[processing.getTdEffective()-1 -i]*=factor;
            spectrum[processing.getTdEffective()-2 -i]*=factor;
        }
    }

    public void sineSquared() {
        double factor;
        double offset = (180.0 - processing.getSsbSineSquared())/180.0;
        // check the loop conditions in the original code...
        for (int i =0 ;
                i< (processing.getTdEffective()-1) && i >= ((processing.getTdEffective()-1)-spectrum.length) ;
                i+=2){
            factor=Math.pow(Math.sin(i/processing.getTdEffective()*Math.PI*offset),2);
            spectrum[processing.getTdEffective()-1 -i]*=factor;
            spectrum[processing.getTdEffective()-2 -i]*=factor;
        }
    }

    /**
     * performs the guassian apodization: W(x,i)=x*exp(-((i*dw)/(2*sigma))^2))
     * where i*dw gives the time coordinate in (s) and sigma is the line broadening with value 0.1 Hz.
     * The cuteNMR implementation is de facto W(x,i)=x*exp(-(i*dw*2*lb)^2))
     * @throws Exception
     */
    public double [] gaussian() throws Exception {
        return gaussian(0.1);
    }

    /**
     * performs the guassian apodization: F(x,i)=x*exp(-((i*dw)/(2*sigma))^2))
     * where i*dw gives the time coordinate in (s) and sigma is the line broadening (lbGauss) in Hz.
     * The cuteNMR implementation is de facto W(x,i)=x*exp(-(i*dw*2*lb)^2))
     * @param lbGauss
     * @throws Exception
     */
    public double [] gaussian(double lbGauss) throws Exception {

        double [] data = new double[processing.getTdEffective()];
        double time;

        if(lbGauss ==0)
            throw new Exception ("line broadening cannot be zero");
        // original expression: (1.0/std::sqrt(2))/par->lbGauss * (1.0/std::sqrt(2))/par->lbGauss;
        double sigmaFactor=0.5*Math.pow(1/lbGauss,2);

        if (acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i++){
                time=i*processing.getDwellTime();
                data[i]=spectrum[i]*Math.exp(-time*time/(sigmaFactor));
            }
        } else {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i+=2){
                time=i*processing.getDwellTime();
                double factor =Math.exp(-time*time/sigmaFactor);
                data[i]=spectrum[i]*factor;
                data[i+1]=spectrum[i+1]*factor;
            }
        }
        return data;
    }

    /**
     * performs the guassian apodization according to Bramer 2001: F(x,i)=x*exp(-(i*dw*lb)^2)
     * where i*dw gives the time coordinate in (s) and the line broadening (lbGauss) in Hz.
     * This also correspondes to one of the alternative ways of defining the standard gaussian distribution,
     * according to wikipedia.
     * @param lbGauss
     * @throws Exception
     */
    public double [] gaussianBramer2001(double lbGauss) throws Exception {
        double [] data = new double[processing.getTdEffective()];
        double time;
        if(lbGauss ==0)
            throw new Exception ("line broadening cannot be zero");
//        double sigmaFactor=0.5*Math.pow(1/lbGauss,2);

        if (acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i++){
                time=i*processing.getDwellTime();
                data[i]=spectrum[i]*Math.exp(-Math.pow(time*lbGauss,2));
            }
        } else {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i+=2){
                time=i*processing.getDwellTime();
                double factor=Math.exp(-Math.pow(time*lbGauss,2));
                data[i]=spectrum[i]*=factor;
                data[i+1]=spectrum[i+1]*factor;
            }
        }
        return data;
    }

    public void firstPoint(){
        firstPoint(1);
    }
    public void firstPoint(double fpCorrection){
        spectrum[0]*=fpCorrection;
        if(!acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            spectrum[1]*=fpCorrection;
        }
    }

    public void zeroFill(){
        zeroFill((int) Math.floor(spectrum.length/5));
    }

    /**
     * This operation adds a determined set of zeros to the fid tail.
     * @param numberZeros
     */
    private void zeroFill(int numberZeros) {
        double[] tmp;
        //TODO Check if the sequencial and simultaneous stuff is right...
        if (acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            tmp= new double[spectrum.length+numberZeros];
        }   else {
            tmp= new double[spectrum.length+numberZeros*2];
        }
        System.arraycopy(spectrum, 0, tmp, 0, spectrum.length);
        for (int i=spectrum.length; i < tmp.length; i++)
            tmp[i]=0;
        spectrum=tmp;
    }
}
