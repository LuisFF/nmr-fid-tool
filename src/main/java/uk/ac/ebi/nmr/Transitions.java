/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.ebi.nmr;

import Jama.Matrix;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author ldpf
 */
public class Transitions {

    private static Double[] transProb;
    private static Integer[] lfs;
    private static Integer[] mfs;
    private static List<Double> frequenciesSpectrum;
    private static List<Double> gintensities;
    private int ik;

    public Transitions(int sizeMax) {
        this.transProb = new Double[sizeMax * sizeMax];
        this.lfs = new Integer[sizeMax * sizeMax];
        this.mfs = new Integer[sizeMax * sizeMax];
        this.frequenciesSpectrum = new ArrayList<Double>();
        this.gintensities = new ArrayList<Double>();

    }
    // this class only changes objects: transProb, 
    // frequenciesSpectrum, 
    // gintensities

    public void computeTransitions(int ifreqSpec,
            int nBasis,
            int nSpins,
            Integer istartPrev,
            int istart,
            Integer sizePrev,
            int size,
            Integer[] basis,
            List<Double> eValues,
            Matrix eVectors,
            List<Double> eValuesPrev,
            Matrix eVectorsPrev) {

        int i = 0;
        int ja, jb;

        int nfs = 0;
//        int sizePrev = binomialCoef.get(istartPrev - 1);
//        int istart = start.get(istartPrev - 1);
        int nDim = sizePrev * size;
        int l, m;
        
        
        for (i = 0; i < nDim; i++) {
            transProb[i] = 0.0;
        }
        /* printf("Begin transition propability calculation\n");*/
        nfs = 0;

        for (l = istartPrev; l < istartPrev + sizePrev; l++) {
            for (m = istart; m < istart + size; m++) {
                int is;
                is = 0;
                for (i = 0; i < nSpins; i++) {
                    is += basis[l + i * nBasis] * basis[m + i * nBasis];
                }
                if (is == nSpins - 2) {
                    lfs[nfs] = l - istartPrev;
                    mfs[nfs] = m - istart;
                    nfs++;
                }
            }
        }
        /* printf("End lfs/mfs calculation nfs %d\n",nfs);*/
//#ifdef ENABLE_OMP
//#pragma omp parallel for private(ja,jb,ik)  
//#endif
        for (ja = 0; ja < sizePrev; ja++) {
            for (jb = 0; jb < size; jb++) {
                double dum = 0.0;
                for (ik = 0; ik < nfs; ik++) {
                    dum += eVectorsPrev.get(lfs[ik], ja)
                            * eVectors.get(mfs[ik], jb);
                }
                transProb[ja * size + jb] = dum * dum;
            }
        }

        /* printf("End transProb calculations\n");*/
        /* transition frequencies*/
        this.ik = ifreqSpec;
        /* printf("Begin gintensities calculations\n");*/
//#ifdef ENABLE_OMP
//#pragma omp parallel for private(ja,jb)  
//#endif
        
        for (ja = 0; ja < sizePrev; ja++) {
            for (jb = 0; jb < size; jb++) {
                int ii = jb + ja * size;
                System.out.println(transProb[ii]+ " "+ Math.abs(eValues.get(jb) - eValuesPrev.get(ja)));
                if (transProb[ii] > .0001) {
                    frequenciesSpectrum.add(Math.abs(eValues.get(jb) - eValuesPrev.get(ja)));
                    gintensities.add(transProb[ii]);
                    ik++;
                }
                    
//                    System.out.println(eValues.get(jb) + " "
//                            + eValues.get(ja) + " "
//                            + frequenciesSpectrum.get(frequenciesSpectrum.size() - 1) + " "
//                            + gintensities.get(gintensities.size() - 1));
//                
            }
        }
        /* printf("End gintensities calculations\n");*/

    }

    public int getIk() {
        return ik;
    }

    public static Integer[] getLfs() {
        return lfs;
    }

    public static Integer[] getMfs() {
        return mfs;
    }

    public static Double[] getTransProb() {
        return transProb;
    }

    public static List<Double> getFrequenciesSpectrum() {
        return frequenciesSpectrum;
    }

    public static List<Double> getGintensities() {
        return gintensities;
    }
}