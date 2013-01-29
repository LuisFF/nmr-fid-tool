package uk.ac.ebi.nmr.fid;

import uk.ac.ebi.nmr.fid.tools.ApodizationTool;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 14/01/2013
 * Time: 14:02
 * To change this template use File | Settings | File Templates.
 */
public class FourierTransformTool {

    private static Acqu acquisition;
    private static Proc processing;
    private static Fid fid;

    private static double[] data;               //proc_buffer   apodization::transform::do_fft






    public FourierTransformTool() {
    }

    public FourierTransformTool(Fid fid, Acqu acquisition) {
        Proc processing = null;
        try {
            processing = new Proc(acquisition);
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        new FourierTransformTool(fid, acquisition, processing);

    }

    public FourierTransformTool(Fid fid, Acqu acquisition, Proc processing) {

        this.fid=fid;
        this.acquisition=acquisition;
        this.processing= processing;

        // instanciating the array where the fourier transformed spectra will be stored....
        switch (acquisition.getAcquisitionMode()) {
            case DISP:
            case SIMULTANIOUS:
                data =new double[2*processing.getTransformSize()];// allocate space for processing
                break;
            case SEQUENTIAL:
                data =new double[4*processing.getTransformSize()];// allocate space for processing
                break;
            default:
                break;
        }

        // perform a left or right shift of the data (ignoring the corresponding portion of the data)
        // the code from cuteNMR was simplified
        for (int i =0; i< (processing.getTdEffective()- Math.abs(processing.getShift())); i++ ) {
                int dataIndex = (processing.getShift()>=0)?
                        i :                             // start in the correct order
                        i - processing.getShift();      // or shift the placement of the data to the right
                int fidIndex = (processing.getShift()>=0)?
                        i + processing.getShift() :     // shift the data to the left (ignore the first indexes)
                        i;                              // start in the correct order
            switch (acquisition.getFidType()){
                case INT32:
                    data[dataIndex]=(double) fid.getData()[fidIndex];
                    break;
                case DOUBLE:
                    data[dataIndex]=fid.getData()[fidIndex];
                    break;
                case FLOAT:
                    break;
                case INT16:
                    break;
                default:
                    break;
            }
        }

        // for sequential data apply Redfield trick: multiply data by 1 -1 -1 1
        if(acquisition.getAcquisitionMode().equals(Acqu.AcquisitionMode.SEQUENTIAL)){
            for (int i = 0; i < processing.getTdEffective(); i+=4) {
                data[i+1] = -(data[i+1]);
                data[i+2] = -(data[i+2]);
            }
        }

        //// try to understand this bit of code!!!!
        //nonsense case if SI is set to be less than TD/2
        int td = (processing.getTdEffective() > 2 * processing.getTransformSize())?
                2 * processing.getTransformSize() :
                processing.getTdEffective();

        // move the data from position i to position 2*i, why?
        if(acquisition.getAcquisitionMode().equals(Acqu.AcquisitionMode.SEQUENTIAL)){
            for (int i =td -1; i>0; i-- ){
                data[2*i] = data[i];
                data[i]=0;
            }
        }
        /////////////////////////////////////////////



    }

    public double [] computeFTT(){
        int signals;
        processing.setLineBroadning(0.3);
        /// applyWindowFunctions //window function need to be applied before FT
        ApodizationTool apodization = new ApodizationTool(data, acquisition.getAcquisitionMode(), processing);
        ApodizationTool apodization2 = new ApodizationTool(data, acquisition.getAcquisitionMode(), processing);
        try {
            apodization.exponential();
            double[] data1 = apodization.getSpectrum();
//            apodization2.gaussian(0.3);
            double[] data2 = apodization2.getSpectrum();
            double[] difference = new double[data1.length];
            for (int i =0 ; i< data1.length; i++){
                difference[i]= data1[i]-data2[i];
            }

            return data1;

        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return apodization.getSpectrum();

        // circular_lsp //if manually chosen, do cls // no used in cuteNMR

//        if (acquisition.getAcquisitionMode().equals(Acqu.AcquisitionMode.SIMULTANIOUS) ||
//                acquisition.getAcquisitionMode().equals(Acqu.AcquisitionMode.DISP)){
//            signals = processing.getTransformSize();
//        }
//
//        if (acquisition.getAcquisitionMode().equals(Acqu.AcquisitionMode.SEQUENTIAL)){
//            signals = 2*processing.getTransformSize();
//        }

//        return null;

    }

    /**
     * translating the fft from the cuteNMR...
     * @return
     */

    public void fft(boolean isForward) {
        int j =1;
        int numberOfPairs=2*processing.getTransformSize();
        int sign = (isForward)? 1 : -1;

        // I have seen this somewhere
        // perhaps I should do a method out of it
        System.out.println(data.length+" "+ numberOfPairs);
        for (int i = 1; i<numberOfPairs; i+=2){
            // inverse the order of the i and j...?
            if( j>i){
                double tmp = data[j-1];
                data[j-1]=data[i-1];
                data[i-1]=tmp;
                tmp=data[j];
                data[j]=data[i];
                data[i]=tmp;
            }
            int m =numberOfPairs/2;
            while (m >= 2 && j > m){
                j=j-m;
                m=m/2;
            }
            j=j+m;

        }
        int mMax=2;
        // changed a bit the calculation of the theta...
        while (numberOfPairs > mMax){
            int iterationStep = 2*mMax;
            double theta= sign*2*Math.PI/mMax;
            double wtmp = Math.sin(0.5*theta);
            double wpr = -2*wtmp*wtmp;
            double wpi = Math.sin(theta);
            double wi=0.0;
            double wr=1.0;
            for(int i=1; i<mMax;i+=2){
                for (int m = i; j<numberOfPairs; i+=iterationStep){
                    j=m+mMax;
                    // position j-1 holds the real point whereas the j holds the imaginary;
                    // thus the wr and wi is a way to weight the corresponding point?
                    double realTmp= wr * data[j-1]-wi*data[j];
                    double imaginaryTmp= wr * data[j]+ wi*data[j-1];
                    data[j-1]=data[m-1]-realTmp;
                    data[j]=data[m]-imaginaryTmp;
                    data[m-1]+=realTmp;
                    data[m]+=imaginaryTmp;
                }
                wtmp=wr;
                wr=wtmp*wr-wi*wpi+wr;
                wi=wi*wpr+wtmp*wi+wi;
            }
            mMax=iterationStep;
        }
        System.out.println("Finished FTT");
    }


    /**
     * The reference bitreverse function.
     */
    private static int bitreverseReference(int j, int nu) {
        int j2;
        int j1 = j;
        int k = 0;
        for (int i = 1; i <= nu; i++) {
            j2 = j1 / 2;
            k = 2 * k + j1 - 2 * j2;
            j1 = j2;
        }
        return k;
    }


    public static Acqu getAcquisition() {
        return acquisition;
    }

    public static Proc getProcessing() {
        return processing;
    }

    public static Fid getFid() {
        return fid;
    }

    public static double[] getData() {
        return data;
    }

    public static void setData(double[] data) {
        FourierTransformTool.data = data;
    }

    public static void setData(int index, double point) {
        FourierTransformTool.data[index] = point;
    }


}
