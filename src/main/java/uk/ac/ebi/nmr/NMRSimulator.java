/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.ebi.nmr;

import Jama.Matrix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 *
 * @author ldpf
 */
public class NMRSimulator {

    

    private int NMAXGROUP = 10;
    private static List<Double> x;
    private static List<Double> y;
    // final variables to be calculated
    private static List<Double> frequencies;
    
    private static List<Double> ppmJCouplings;
    // basis stuff
    private static Integer[] basis;
    private static List<Integer> start;
    private static List<Integer> binomialCoef;
    // counting entities
    private static int nFrequencies;
    private static int nCoup;
    private static int nSpins;
    private static int nBasis;

    public NMRSimulator() {
    }

    public void computeNMRSpectrum(double operatingFrequency,
            int nGroups,
            List<Integer> numberOfSpins,
            List<Double> chemichalShifts,
            Double[][] jCouplings,
            int level) throws Exception {

        if (jCouplings == null) {
            computeNMRSpectrumNonCoupled(nGroups, numberOfSpins, chemichalShifts);
        } else if (level == -1) {
            computeNMRSpectrumOldMethod(operatingFrequency, nGroups, numberOfSpins, chemichalShifts, jCouplings);
        } else if (level == 2) {
            computeNMRSpectrumByBlock(operatingFrequency, nGroups, numberOfSpins, chemichalShifts, jCouplings);
        } else if (level == 1) {
            computeNMRSpectrumWeaklyCoupled(operatingFrequency, nGroups, numberOfSpins, chemichalShifts, jCouplings);
        } else {
            computeNMRSpectrumNonCoupled(nGroups, numberOfSpins, chemichalShifts);
        }

    }

    private void computeNMRSpectrumNonCoupled(int nGroups, List<Integer> numberOfSpins, List<Double> chemichalShifts) throws Exception {
        this.nSpins = 0;

        x= new ArrayList<Double>();
        y = new ArrayList<Double>();
        nFrequencies = 0;

        setNumberSpins(numberOfSpins);



        if (nSpins < 1) {
            throw new Exception("Not enough number of spins");
        }

        this.nFrequencies = nSpins;


        // expand the spin systems to single atoms attributing the same shift
        // and intensity to atoms that belong to the same spin system
        for (int i = 0; i < chemichalShifts.size(); i++) {
            for (int j = 0; j < numberOfSpins.get(i); j++) {
                x.add(chemichalShifts.get(i));
                y.add(1.0);
            }
        }
    }

    private void setNumberSpins(List<Integer> numberOfSpins) {
        // the following method was simplified 
//        getNumberOfSpinsCouples(nGroups, numberOfSpins, nSpins, nCoup);
        // count the number of spins
        for (Integer n : numberOfSpins) {
            nSpins += n;
        }
        // calculate all possible spin combinations
        // in the original source this cas calculated using for loops
        // I think that the number of couplings possible is 2^nSpins
        this.nCoup = 1 << nSpins;

        /////////
    }

    private void computeNMRSpectrumOldMethod(double operatingFrequency, int nGroups, List<Integer> numberOfSpins, List<Double> chemichalShifts, Double[][] jCouplings) {
        throw new UnsupportedOperationException("Not yet implemented");
    }

    private void computeNMRSpectrumByBlock(double operatingFrequency, int nGroups, List<Integer> numberOfSpins, List<Double> chemichalShifts, Double[][] jCouplings) throws Exception {
        nSpins = 0;
        int i = 0;
        int nDim = 0;
        nCoup = 0;


//        frequenciesSpectrum = new ArrayList<Double>();
//        gintensities = new ArrayList<Double>();


        Double[] transProb;
        Double[] hamilt;
        Matrix eVectors;
        List<Double> eValues;
        Matrix eVectorsPrev;
        List<Double> eValuesPrev;
        List<Integer> lfs;
        List<Integer> mfs;
        nFrequencies = 0;

        int sizeMax;
        int sizeTrans;
        int ifreqSpec = 0;

//	*n = 0;
//	*X = NULL;
//	*Y = NULL;

        setNumberSpins(numberOfSpins);

        /* printf("# Spins = %d\n",nSpins);*/
        if (nSpins < 1) {
            throw new Exception("Nof enough spin numbers");
        }
        if (nSpins > NMAXGROUP + 2) {
            throw new Exception("The maximum number of spins is " + (NMAXGROUP + 2));
        }

        frequencies = new ArrayList<Double>(nSpins);
        ppmJCouplings = new ArrayList<Double>();

        /* frequency array and coupling matrix*/
        // build an expanded version of the shifts and coupling arrays
        buildFreqAndJ(operatingFrequency, nGroups, numberOfSpins, chemichalShifts, jCouplings);
        /* printf("End build FreqJ\n");*/
        /* basis set */
        // the following operation corresponds to 2^nspins
        // this defines the size of the basis that one has to calculate
        nBasis = 1 << nSpins;

        /* printf("nBasis =%d\n",nBasis);*/
        basis = new Integer[nBasis*nSpins];
        start = new ArrayList<Integer>();
        binomialCoef = new ArrayList<Integer>();

        // initial table?
        initTables();
        // get maximal binomial coef
//        sizeMax = getSizeMax(nSpins, binomialCoef);
        sizeMax = Collections.max(binomialCoef);
        // sizeTrans? sum of the product between the binomial coeff of position i-1 and i
        // this is used to initialize the variables to store the x and y coordinates
        sizeTrans = getSizeTrans();
        System.out.println("nBasis=" + nBasis + " sizeTrans=" + sizeTrans);

        /* printf("nBasis  = %d sizeTrans = %d\n",nBasis,sizeTrans);*/
        /* printf("End initTables\n");*/

        nDim = sizeMax * (sizeMax + 1) / 2;
        //hamilt = new Double[nDim];
        transProb = new Double[sizeMax * sizeMax];// transition probability?
        /* printf("End transProb alloc\n");*/

        eValues = new ArrayList<Double>(sizeMax);
//        eVectors = new ArrayList<List<Double>>();
        // initialize eVectors? Do I need that in Java?
//        if (eVectors) {
//            for (i = 0; i < sizeMax; i++) {
//                eVectors[i] = malloc(sizeMax * sizeof(gdouble));
//            }
//        }
        eValuesPrev = new ArrayList<Double>();
        //eVectorsPrev = new ArrayList<List<Double>>();
        eVectorsPrev = null;
        // again initializing also the list inside the list. Do I need that?
//        if (eVectorsPrev) {
//            for (i = 0; i < sizeMax; i++) {
//                eVectorsPrev[i] = malloc(sizeMax * sizeof(gdouble));
//            }
//        }

        /* printf("Begin alloc freq and gint\n");*/
//        frequenciesSpectrum = new ArrayList<Double>(sizeTrans);
//        gintensities = new ArrayList<Double>(sizeTrans);
        //check if there is enough memory allocated for the intensities...
        // lets give it a shot whith this first
//        if (!gintensities) {
//            gchar tmp[
//            BSIZE
//            ];
//		sprintf(tmp, _("Sorry\n Not enough memory : nBasis  = %d sizeTrans = %d\n"), nBasis, sizeTrans);
//            Message(tmp, _("Warning"), TRUE);
//            return;
//        }
        lfs = new ArrayList<Integer>();
        mfs = new ArrayList<Integer>();
        ifreqSpec = 0;
        Transitions trs = new Transitions(sizeMax);

        for (i = 0; i <= nSpins; i++) {
            int istart = start.get(i);
            int size = binomialCoef.get(i);
            
//            
            ///check on HAMILTONIAN!!! What does it do and what is it good for!!!
            //Hamiltonian is the operator corresponding to the total energy of the system
            hamilt = buildHamiltonianOneBlock(istart, size, nBasis, nDim);
            /* prgintMatInf(0, size, hamilt, "H");*/
            /*diagonaliseJacobiOneBlock(size, hamilt, eValues, eVectors);*/

            /* printf("Diag size = %d\n",size);*/
            // calculate the eigen values and eigen vectors of the hamiltononian
            // I have to build a matrix out of the hamiltonian first, using the size
            
            
            Matrix hamiltM = getHamiltMatrix(hamilt, size);
            // the eifen also scales the hamilt 
            // That's probably the reason why I get things a bit messed up...
//            eigen(hamilt, size, eValues, eVectors);
            // new way to get the eigenvalues and the eigen vectors
            // let's see if this works
            Matrix eigen = hamiltM.eig().getD();
            // In principle the eigen values should be the diagonal of the matrix
            // but you never know...
            // So I have two alternatives. 
//            for (int j = 0; j < eigen.getRowDimension(); j++) {
//                for (int k = 0; k < eigen.getColumnDimension(); k++) {
//                    if(Math.abs(eigen.get(j,k))>0.001)
//                        eValues.add(eigen.get(j,k));
//                }
//            }
            for (int j = 0; j < eigen.getRowDimension(); j++) {
                        eValues.add(eigen.get(j,j));
            }
            
            eVectors = hamiltM.eig().getV();
            System.out.println(eValues.size()+ " "+ eVectors.rank());
            // previous approach by Allouche
//            QL eigen = new QL(Arrays.asList(hamilt),size,sizeMax);
//            eigen.computeEigen();
//            eValues = Arrays.asList(eigen.geteVals());
//            eVectors = new Matrix(eigen.getvMatrix());

            
            
            if (i != 0) {
                /* printf("Begin transition calculations\n");*/
                
                trs.computeTransitions(ifreqSpec, nBasis,nSpins,start.get(i-1), istart, binomialCoef.get(i-1), size,
					basis,eValues, eVectors,eValuesPrev,  
					eVectorsPrev);
                ifreqSpec = trs.getIk();
                
                
//					eValuesPrev,  
//					eVectorsPrev,
//					, frequenciesSpectrum, gintensities,transProb,lfs,mfs);
            }
//            /* swith prev and current eigenvalues and eigenvectors */
            List<Double> e = eValuesPrev;
            eValuesPrev = eValues;
            eValues = e;
            
            Matrix v = eVectorsPrev;
            eVectorsPrev = eVectors;
            eVectors = v;
        }
//        frequencies=null;
//        ppmJCouplings=null;
        hamilt=null;
        eValues=null;
        eValuesPrev=null;
        eVectors=null;
        eVectorsPrev=null;
        transProb=null;
        
//        if (frequencies) {
//            free(frequencies);
//        }
//        if (ppmJCouplings) {
//            free(ppmJCouplings);
//        }
//        if (hamilt) {
//            free(hamilt);
//        }
//        if (lfs) {
//            free(lfs);
//        }
//        if (mfs) {
//            free(mfs);
//        }
//        /* prgint Hamiltonian */
//
//        if (eVectorsPrev) {
//            for (i = 0; i < sizeMax; i++) {
//                if (eVectorsPrev[i]) {
//                    free(eVectorsPrev[i]);
//                }
//            }
//            free(eVectorsPrev);
//        }
//        if (eVectors) {
//            for (i = 0; i < sizeMax; i++) {
//                if (eVectors[i]) {
//                    free(eVectors[i]);
//                }
//            }
//            free(eVectors);
//        }
//        if (eValues) {
//            free(eValues);
//        }
//        if (eValuesPrev) {
//            free(eValuesPrev);
//        }
//        if (transProb) {
//            free(transProb);
//        }
//
        nFrequencies = ifreqSpec;
//        if (nFrequencies > 0) {
//            /*
//            gint i;
//            for(i=0;i<nFrequencies;i++) 
//            {
//            printf("%f %f\n",frequenciesSpectrum[i],gintensities[i]);
//            }
//             */
//            /*
//            nFrequencies =  removeIdenticalFrequencies(nFrequencies, frequenciesSpectrum, gintensities);
//             */
//            if (nFrequencies < 1) {
//                frequenciesSpectrum=null;
//                gintensities=null;
////                if (frequenciesSpectrum) {
////                    free(frequenciesSpectrum);
////                }
////                frequenciesSpectrum = NULL;
////                if (gintensities) {
////                    free(gintensities);
////                }
////                gintensities = NULL;
//            } else {
//                
//                frequenciesSpectrum = realloc(frequenciesSpectrum, nFrequencies * sizeof(gdouble));
//                gintensities = realloc(gintensities, nFrequencies * sizeof(gdouble));
//            }
//        }
//
//         * n = nFrequencies;
//         * X = frequenciesSpectrum;
//         * Y = gintensities;
        x = trs.getFrequenciesSpectrum();
        y = trs.getGintensities();
        System.out.println("Hurray!!!");

    }

    private void computeNMRSpectrumWeaklyCoupled(double operatingFrequency,
            int nGroups,
            List<Integer> numberOfSpins,
            List<Double> chemichalShifts,
            Double[][] jCouplings) {
        throw new UnsupportedOperationException("Not yet implemented");
    }


    
    // This seems to be working fine. 
    // I have tested with spin systems bigger than 1
    private void buildFreqAndJ(double operatingFrequency, 
            int nGroups, 
            List<Integer> numberOfSpins, 
            List<Double> chemichalShifts, 
            Double[][] jCouplings) {
        
        for (int i = 0; i < chemichalShifts.size(); i++) {
            // expand the sping systems to single atoms
            // create the same shifts for all the atoms in the current spin system
            // create also the same coupling to other spin systems for the atoms in the current spin systems
//            System.out.println(chemichalShifts.size() + " " + numberOfSpins.size() + " " + i);
            for (int n = 0; n < numberOfSpins.get(i); n++) {
                frequencies.add(chemichalShifts.get(i));
                // define the same couplings for all the atoms in the same group
                // process information columnwise (starting for the 1st atom)
                for (int j = 0; j < i; j++) {
                    //add the coupling between atom j and each atom in the previous spin systems
                    //the previous spin systems also need to be decoupled in single atoms
                    for (int m = 0; m < numberOfSpins.get(j); m++) {
                        if (jCouplings[i][j] != null) {
                            ppmJCouplings.add(jCouplings[i][j] / operatingFrequency);
                        }
                    }
                }
                //put zero to all the atoms in the same spin system
                for (int l = 0; l <= n; l++) {
                    ppmJCouplings.add(0.0);
                }
            }
        }
    }

    private void initTables() {
        /* binomialCoef binomial coefficients. nSpins+1 elements */
        /* bIndex element i is locn of fn i in basis[] */
        /* start indexes to successive Fz blocks in basis[]. nSpins+1 elements */
        /* basis : coded basis functions  : 2**nSpins*nSpins*/
        int m, k, i, j, ik;
        
        List<Integer> sumSpin;

        binomialCoef.add(1);
        start.add(0);

        //generate a bionomial series (Pascal's triangle)
        // this is to be used for the peak multiplicities stuff
        for (i = 1, k = nSpins, m = 1; i <= nSpins; i++, k--, m++) {
            Integer lastBionomial = binomialCoef.get(binomialCoef.size()-1);

            // check if this generates an integer in java!!!
            binomialCoef.add((lastBionomial * k) / m);
//		binomialCoef[i] = (binomialCoef[i-1] * k) / m;

            // sum the last start with the last binomial???
            // 0+1; 1+2; 3+1;
            start.add(start.get(start.size()-1) + lastBionomial);
//		start[i] = start[i-1] + binomialCoef[i-1];
        }
        // Bitwise left shift (I guess it's the same thing in Java)
        // 2^nSpins
        // it has been made global variable
        //nBasis = 1 << nSpins;
        /* basis = -1 (alpha) or +1(beta) */
        // The basis is indeed a basis like in EFMS
        // each vector is independet from each other
        // this represents all the possible combinations between atom spins(?)
        for (i = 0; i < nBasis; i++) {
            int a = i;
            //?????
            for (j = 0; j < nSpins; j++) {
                // nSpins -current spin -1
                // spins left -1?
                k = nSpins - j - 1;
                // 2^(spins left -1)
                int b = 1 << (k);
                ik = i + k * nBasis;
                basis[ik]= a / b;
                a = a - b * basis[ik];
                //flag entries that are zero with a -1
                if (basis[ik] == 0) {
                    basis[ik]= -1;
                }
            }
        }
        // The spumSpin is used to define the intensity of each peak combination(?)
        sumSpin = new ArrayList<Integer>();
        //define the intensity of the peaks??????
        for (i = 0; i < nBasis; i++) {
            int is = 0;
            for (j = 0; j < nSpins; j++) {
                is += basis[i + j * nBasis];
            }
            sumSpin.add(is);
        }

        for (i = 0; i < nBasis - 1; i++) {
            k = i;
            // compare the intensities???
            for (j = i + 1; j < nBasis; j++) {
                if (sumSpin.get(k) > sumSpin.get(j)) {
                    k = j;
                }
            }
            // if sumSpin.get(k)>sumSpin.get(j)
            if (k != i) {
                // flip the values of position k with position i in sumSpin
                int t = sumSpin.get(k);
                sumSpin.set(k, sumSpin.get(i));
                sumSpin.set(i, t);
                /////////
                for (j = 0; j < nSpins; j++) {
                    // flip also the values of the basis for all the spins in between???
                    t = basis[i + j * nBasis];
                    basis[i + j * nBasis]= basis[k + j * nBasis];
                    basis[k + j * nBasis]=t;
                }
            }
        }
        //delete the sumSpin variable
        //free(sumSpin);
        // Java way of doing stuff
        sumSpin = null;

    }

    private int getSizeTrans() {
        int size = 0;
        for (int i = 1; i <= nSpins; i++) {
            size += binomialCoef.get(i) * binomialCoef.get(i - 1);
        }
        return size;
    }

    private Double[] buildHamiltonianOneBlock(int istart, int size, int nBasis, int nDim) {
        int i = 0;
        int j = 0;
        int k = 0;
        int l = 0;
        int ja, jb;
        int ii, ik, is;
        double dum;
        int dumja, dumjb;
        Double[] hamilt = new Double[nDim];

        /* diagonal elements of hamiltonian*/
        ii = 0;
        for (i = istart; i < istart + size; i++) {
            dum = 0.0;
            ik = -1;
            for (l = 0; l < nSpins; l++) {
                ik++;
                dum += frequencies.get(l) * basis[i + l * nBasis] / 2;
                for (k = 0; k < l; k++) {
                    dum = dum + ppmJCouplings.get(ik) * basis[i + l * nBasis] * basis[i + k * nBasis] / 4.;
                    ik++;
                }
            }
            ii = ii + i + 1 - istart;
            hamilt[ii - 1]= dum;
            /* printf("hamilt[%d] = %f\n",ii-1, hamilt[ii-1]);*/
        }
        /* printf("End calculate diagonal elements of hamiltonian\n");*/

        /*	calculate off diagonal elements of hamiltonian*/
        ii = -1;
        for (ja = istart + 1; ja < istart + size; ja++) {
            ii = ii + 1;
            for (jb = istart; jb < ja; jb++) {
                ii = ii + 1;
                is = 0;
                dumja = 0;
                dumjb = 0;
                /* to delete ? */
                for (i = 0; i < nSpins; i++) {
                    is += basis[ja + i * nBasis] * basis[jb + i * nBasis];
                    dumja += basis[ja + i * nBasis];
                    dumjb += basis[jb + i * nBasis];
                }
                if ((dumja != dumjb) || (is != nSpins - 4)) {
                    hamilt[ii]= 0.0;
                    /* printf("ii = %d hamilt[ii] = %f\n",ii, hamilt[ii]);*/
                    continue;
                }
                for (i = 0; i < nSpins - 1; i++) {
                    if (basis[ja + i * nBasis] != basis[jb + i * nBasis]) {
                        break;
                    }
                }
                for (j = i + 1; j < nSpins; j++) {
                    if (basis[ja + j * nBasis] != basis[jb + j * nBasis]) {
                        break;
                    }
                }
                hamilt[ii]= 0.5 * ppmJCouplings.get((j + 1) * j / 2 + i);
                /* printf("ii = %d hamilt[ii] = %f\n",ii, hamilt[ii]);*/
            }
        }
        return hamilt;
    }

    private Matrix getHamiltMatrix(Double[] hamilt, int n) {
        Matrix h = new Matrix(n, n);
        int ii = -1;
	for(int i=0;i<n;i++)
	for(int j=0;j<=i;j++)
	{
		ii++;
		h.set(i, j, hamilt[ii]);
	}
	for(int i=0;i<n;i++)
  	for(int j=i+1;j<n;j++)
            h.set(i, j, h.get(j, i));
        return h;
    }

    public static List<Double> getX() {
        return x;
    }

    public static List<Double> getY() {
        return y;
    }

    public static int getnFrequencies() {
        return nFrequencies;
    }
    

    
}
