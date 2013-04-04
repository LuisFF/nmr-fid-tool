package uk.ac.ebi.nmr.fid;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 14/01/2013
 * Time: 14:02
 * To change this template use File | Settings | File Templates.
 */
public class LdpfsFourierTransformTool extends AbstractFastFourierTransformTool implements FastFourierTransformTool {

    public LdpfsFourierTransformTool(Fid fid, Acqu acquisition) throws Exception{                
        this.processing = new Proc(acquisition);        
        this.fid = fid;
    }

    public LdpfsFourierTransformTool(Fid fid, Acqu acquisition, Proc processing) {
        this.fid=fid;
        this.acquisition=acquisition;
        this.processing=processing;
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

    /**
     * This is where the implementation of each fft package go.
     * @param apodizedData
     * @return 
     */
    double[] implementedFFT(double[] apodizedData) {
        throw new UnsupportedOperationException("Not yet implemented");
    }


}
