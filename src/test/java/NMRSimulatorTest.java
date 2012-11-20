/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import junit.framework.TestCase;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import uk.ac.ebi.nmr.NMRSimulator;
import uk.ac.ebi.nmr.SpectrumGenerator;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author ldpf
 */
public class NMRSimulatorTest extends TestCase {

    NMRSimulator nmrS;

    /**
     * Create the test case
     *
     */
//    public NMRSimulatorTest( String testName )
//    {
//        super( testName );
//        this.nmrS = new NMRSimulator();
//    }
    public NMRSimulatorTest() {
        this.nmrS = new NMRSimulator();
    }
//
//    /**
//     * @return the suite of tests being tested
//     */
//    public static Test suite()
//    {
//        return new TestSuite( NMRSimulator.class );
//    }

    /**
     * Rigourous Test :-)
     */
    public void testNMRSimulator() {
        double operatingFrequency = 300;
//        double operatingFrequency = 500;
//        int nGroups = 6;
        int nGroups = 2;
//        int nGroups = 5;
        List<Integer> numberOfSpins = new ArrayList<Integer>();
        List<Double> chemichalShifts = new ArrayList<Double>();
        Double[][] jCouplings = new Double[nGroups][nGroups];
        double[] x;
        double[] y;
        double[][] xypot;
        int n;
        int i;

        // TEST in the cpp src 

//        numberOfSpins.add(1);
//        numberOfSpins.add(1);
//        numberOfSpins.add(1);
//        numberOfSpins.add(1);
//        numberOfSpins.add(1);
//        numberOfSpins.add(1);
//        
//
//
//        chemichalShifts.add(2.14);
//        chemichalShifts.add(2.13);
//        chemichalShifts.add(1.68);
//        chemichalShifts.add(1.68);
//        chemichalShifts.add(2.63);
//        chemichalShifts.add(2.63);



//	jCouplings.add(Arrays.asList(Double.valueOf(0)));
//	jCouplings.add(Arrays.asList(10.23,10.245));
//	jCouplings.add(Arrays.asList(10.23,10.245,0.0));
//	jCouplings.add(Arrays.asList(0.0,0.0,9.52,9.58));
//	jCouplings.add(Arrays.asList(0.0,0.0,9.5,9.5,0.0));

//        jCouplings[1][0] = 0.1;
//
//        jCouplings[2][0] = 10.23;
//        jCouplings[2][1] = 10.245;
//
//        jCouplings[3][0] = 10.23;
//        jCouplings[3][1] = 10.245;
//        jCouplings[3][2] = 0.0;
//
//        jCouplings[4][0] = 0.0;
//        jCouplings[4][1] = 0.0;
//        jCouplings[4][2] = 9.52;
//        jCouplings[4][3] = 9.58;
//
//        jCouplings[5][0] = 0.0;
//        jCouplings[5][1] = 0.0;
//        jCouplings[5][2] = 9.5;
//        jCouplings[5][3] = 9.5;
//        jCouplings[5][4] = 0.0;

        // TEST default spectra in gabedit
        numberOfSpins.add(1);
        numberOfSpins.add(2);
        
        chemichalShifts.add(1.5);
        chemichalShifts.add(3.0);
        
        jCouplings[1][0] = 10.0;

        // TEST acetovanillone
//        numberOfSpins.add(1);
//        numberOfSpins.add(1);
//        numberOfSpins.add(1);
//        numberOfSpins.add(3);
//        numberOfSpins.add(3);
//        
//        
//
//
//        chemichalShifts.add(8.19985);
//        chemichalShifts.add(8.66735);
//        chemichalShifts.add(7.13665);
//        chemichalShifts.add(2.90988);
//        chemichalShifts.add(4.53472);
//        
//        
//        jCouplings[1][0] = 1.49;
//
//        jCouplings[2][0] = -0.279;
//        jCouplings[2][1] = 9.37;
//
//        jCouplings[3][0] = 0.375;
//        jCouplings[3][1] = -0.216;
//        jCouplings[3][2] = -0.298;
//
//        jCouplings[4][0] = -0.0961;
//        jCouplings[4][1] = -0.271;
//        jCouplings[4][2] = -0.258;
//        jCouplings[4][3] = -0.0843;

        // TEST alanine
//        numberOfSpins.add(1);
//        numberOfSpins.add(2);
//
//        chemichalShifts.add(1.5);
//        chemichalShifts.add(3.0);
//
//        jCouplings[1][0] = 10.0;

        try {
//            nmrS.computeNMRSpectrum(operatingFrequency, nGroups, numberOfSpins, chemichalShifts, jCouplings, -1);
            nmrS.computeNMRSpectrum(operatingFrequency, nGroups, numberOfSpins, chemichalShifts, jCouplings, 2);

            int nbOfStates = nmrS.getnFrequencies();
            List<Double> energies = nmrS.getX();
            List<Double> intensities = nmrS.getY();
            System.out.println(energies.size() + " " + intensities.size());

            SpectrumGenerator spectG = new SpectrumGenerator(energies, intensities);
            spectG.computeCurvePoints();
            spectG.computePeaks();
            spectG.setYmax2One();




            XYSeries series = new XYSeries("NMR spectra");
            int count = 0;
            for (int j = 0; j < spectG.getX_curve().length; j++) {
                series.add(spectG.getX_curve()[j], spectG.getY_curve()[j]);
                count++;
                if (count > 10) {
                    System.out.println(spectG.getX_curve()[j] + " " + spectG.getY_curve()[j]);
                    count = 0;
                }
            }
            count = 0;

            for (int j = 0; j < spectG.getX_peaks().length; j++) {
                series.add(spectG.getX_peaks()[j], spectG.getY_peaks()[j]);
                count++;
                if (count > 10) {
                    System.out.println(spectG.getX_peaks()[j] + " " + spectG.getY_peaks()[j]);
                    count = 0;
                }
            }
            System.out.println("Number of points: " + spectG.getX_peaks().length + " "
                    + spectG.getY_peaks().length);


            //         Add the series to your data set
            XYSeriesCollection dataset = new XYSeriesCollection();
            dataset.addSeries(series);
            //         Generate the graph
            JFreeChart chart = ChartFactory.createXYLineChart("1H NMR", // Title
                    "x-axis", // x-axis Label
                    "y-axis", // y-axis Label
                    dataset, // Dataset
                    PlotOrientation.VERTICAL, // Plot Orientation
                    true, // Show Legend
                    true, // Use tooltips
                    false // Configure chart to generate URLs?
                    );
            try {
                ChartUtilities.saveChartAsJPEG(new File("/Users/ldpf/Downloads/Chart.jpg"), chart, 1500,
                        1300);
            } catch (IOException e) {
                System.err.println("Problem occurred creating chart.");
            }



        } catch (Exception ex) {
            Logger.getLogger(NMRSimulatorTest.class.getName()).log(Level.SEVERE, null, ex);
        }

        assertTrue(true);
    }
}
