/*
 * Copyright (c) 2013. EMBL, European Bioinformatics Institute
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


package uk.ac.ebi.nmr.fid.tools.apodization;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;
import uk.ac.ebi.nmr.fid.Spectrum;
import uk.ac.ebi.nmr.fid.io.FidReader;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;

/**
 * Tests for various window functions
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 22/01/2013
 * Time: 16:38
 *
 */

public class ApodizatorTest {

    static Acqu acquisition;
    static Proc processing;
    static Spectrum spectrum;

    static double [] fid;
    static double [] fidApodizedEM;
    static double [] fidApodizedGM;
    static double [] fidApodizedLGM;
    static double [] fidApodizedTRAF;
    static double [] fidApodizedTRAFS;
    static double [] realFTObserved;
    static double [] imaginaryFTObserved;

    XYSeries data;
    static XYSeries dataR;

    @BeforeClass
    public static void loadExternalData () {


        try {
            // the dsp in the name corresponds to the aquisition mode and to the the group delay correction
            fid= (double[]) loadSerializedObject(FidReader.class.getClassLoader()
                    .getResourceAsStream("data/simulated/hypothetical-compound/1h/fid-sim-raw-dsp.ser"));
            fidApodizedGM= (double[]) loadSerializedObject(FidReader.class.getClassLoader()
                    .getResourceAsStream("data/simulated/hypothetical-compound/1h/fid-apgm-sim-dsp.ser"));
            fidApodizedEM= (double[]) loadSerializedObject(FidReader.class.getClassLoader()
                    .getResourceAsStream("data/simulated/hypothetical-compound/1h/fid-apem-sim-dsp.ser"));
            fidApodizedLGM= (double[]) loadSerializedObject(FidReader.class.getClassLoader()
                    .getResourceAsStream("data/simulated/hypothetical-compound/1h/fid-aplgm-sim-dsp.ser"));
            fidApodizedTRAF= (double[]) loadSerializedObject(FidReader.class.getClassLoader()
                    .getResourceAsStream("data/simulated/hypothetical-compound/1h/fid-aptrafm-sim-dsp.ser"));
            fidApodizedTRAFS= (double[]) loadSerializedObject(FidReader.class.getClassLoader()
                    .getResourceAsStream("data/simulated/hypothetical-compound/1h/fid-aptrafsm-sim-dsp.ser"));

            acquisition = new Acqu(Acqu.Spectrometer.BRUKER);
            // The simulated data is for DISP
            acquisition.setAcquisitionMode(Acqu.AcquisitionMode.DISP);
            System.out.println(acquisition.getAcquisitionMode().toString());
            acquisition.setAquiredPoints(fid.length);
//            acquisition.setTransmiterFreq(8 + 1 / 3); // for bmse000109 is app. 500
//            acquisition.setSpectralWidth(6.25); // for bmse000109 is app. 13
            acquisition.setTransmiterFreq(1); // for bmse000109 is app. 500
            acquisition.setSpectralWidth(4000); // for bmse000109 is app. 13
            acquisition.setDspFirmware(10);
            acquisition.setDspDecimation(1);// pass through the dsp phase correction without changing the data
//            processing = new Proc(acquisition);
//            processing.setLineBroadening(0.3);
//            processing.setGbFactor(0.03);
//            System.out.println("DW: " +processing.getDwellTime());
//            System.out.println("LB: " +processing.getLineBroadening());
            spectrum = new Spectrum(fid,acquisition);
            dataR = new XYSeries("original fid");
            for (int i = 0 ; i< fid.length; i++)
                dataR.add(i,fid[i]);

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (ClassNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


    }

    private static Object loadSerializedObject(InputStream inputStream) throws IOException, ClassNotFoundException {
        long startTime = System.currentTimeMillis();
        ObjectInputStream ois = new ObjectInputStream(inputStream);
        Object serializedOject = ois.readObject();
        long endTime   = System.currentTimeMillis();
        System.out.println("Time reading object in "+inputStream.toString()+":\t\t"+(endTime - startTime)+" ms");
        return serializedOject;
    }

    @Test

    public void TestExponentialApodizatorDISP(){

        Assert.assertNotNull("FID was not properly read",spectrum.getFid());

        data = new XYSeries("data");
        Apodizator expApodizator= new ExponentialApodizator(spectrum);

        try {
            Spectrum apodizedSpectrum = expApodizator.calculate(0.3);
            boolean controlSTERR = false;
            Assert.assertArrayEquals("Numerical deviation above 1E-12", fidApodizedEM, apodizedSpectrum.getFid(),1E-12);
            // consider removing this bit below

            for (int i =0; i< apodizedSpectrum.getFid().length;i++){
                data.add(i, apodizedSpectrum.getFid()[i]);
                if(!(controlSTERR) && Math.abs(apodizedSpectrum.getFid()[i]-fidApodizedEM[i])>1E-12) {
                    System.err.println("Exponential apodization with deviation larger than 1E-12");
                    System.err.println(i + " " + fidApodizedEM[i] + " " + apodizedSpectrum.getFid()[i]);
                    controlSTERR=true;
                }
            }
            } catch (Exception e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            XYSeriesCollection dataset = new XYSeriesCollection(dataR);
            dataset.addSeries(data);
            JFreeChart chart = createChart(dataset);
            try {
                ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart-exp.png"), chart, 864 ,1152);
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
    }

    @Test
    public void TestGaussianApodizatorDISP(){
        data = new XYSeries("data");

        Apodizator gaussApodizator= new GaussianApodizator(spectrum);

        try {
            Spectrum apodizedSpectrum = gaussApodizator.calculate(0.3);
            Assert.assertArrayEquals("Numerical deviation above 1E-12", fidApodizedGM, apodizedSpectrum.getFid(),1E-12);

            // consider removing this bit below
            boolean controlSTERR = false;
            for (int i =0; i< apodizedSpectrum.getFid().length;i++){
                data.add(i,apodizedSpectrum.getFid()[i]);
                if(Math.abs(apodizedSpectrum.getFid()[i]-fidApodizedGM[i])>1E-12) {
                    if(!controlSTERR){
                        System.err.println("Gaussian apodization with deviation larger than 1E-12");
                        System.err.println(i+" "+ fidApodizedGM[i]+" "+ apodizedSpectrum.getFid()[i]);
                        controlSTERR=true;
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        XYSeriesCollection dataset = new XYSeriesCollection(dataR);
        dataset.addSeries(data);
        JFreeChart chart = createChart(dataset);
        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart-gauss.png"), chart, 864 ,1152);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    @Test
    public void TestLorentzGausApodizatorDISP(){
        data = new XYSeries("data");
        Apodizator lgApodizator= new LorentzGausApodizatior(spectrum);
        spectrum.getProc().setGbFactor(0.03);
        //changes in the definition of the TdEffective will affect LotentzGaus, TRAF and TRAFS Apodizators
        spectrum.getProc().setTdEffective(acquisition.getAquiredPoints()/2);

        try {
            Spectrum apodizedSpectrum = lgApodizator.calculate(-0.3);
            Assert.assertArrayEquals("Numerical deviation above 1E-4", fidApodizedLGM, apodizedSpectrum.getFid(),1E-4);
            // consider removing this bit below
            boolean controlSTERR = false;
            int count=0;
            for (int i =0; i< apodizedSpectrum.getFid().length;i++){
                data.add(i,apodizedSpectrum.getFid()[i]);
                if(Math.abs(apodizedSpectrum.getFid()[i]-fidApodizedLGM[i])>1E-5) {
                    if(!(controlSTERR)){
                    System.err.println("Lorentz-Gauss apodization with deviation larger than 1E-5");
                    System.err.println(i+" "+fidApodizedLGM[i]+" "+ apodizedSpectrum.getFid()[i]);
                        controlSTERR=true;
                    }
                    count++;
                }
            }
            System.out.println("Lorentz-Gauss apodization deviations above 1E-5: "+count);
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        XYSeriesCollection dataset = new XYSeriesCollection(dataR);
        dataset.addSeries(data);
        JFreeChart chart = createChart(dataset);
        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart-lgauss.png"), chart, 864 ,1152);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    @Test
    public void TestTRAFApodizatorSequential(){
        data = new XYSeries("data");
        Apodizator tarfApodizator = new TrafApodizator(spectrum);

        try {
            Spectrum apodizedSpectrum = tarfApodizator.calculate(0.06);
            Assert.assertArrayEquals("Numerical deviation above 1E-5", fidApodizedTRAF, apodizedSpectrum.getFid(),1E-5);
            // consider removing this bit below
            boolean controlSTERR = false;
            int count = 0;
            for (int i =0; i< apodizedSpectrum.getFid().length;i++){
                data.add(i,apodizedSpectrum.getFid()[i]);
                if(Math.abs(apodizedSpectrum.getFid()[i]-fidApodizedTRAF[i])>1E-6) {
                    if(!(controlSTERR)){
                    System.err.println("TRAF apodization with deviation larger than 1E-6");
                    System.err.println(i+" "+fidApodizedTRAF[i]+" "+ apodizedSpectrum.getFid()[i]);
                    controlSTERR=true;
                    }
                    count++;
                }
            }
            System.out.println("TRAF apodization deviations above 1E-6: "+count);
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        XYSeriesCollection dataset = new XYSeriesCollection(dataR);

        dataset.addSeries(data);
        JFreeChart chart = createChart(dataset);
        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart-traf.png"), chart, 864 ,1152);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    @Test
    public void TestTRAFSApodizatorDISP(){
        data = new XYSeries("data");
        Apodizator tarfsApodizator= new TrafsApodizator(spectrum);

        try {
            Spectrum apodizedSpectrum = tarfsApodizator.calculate(0.06);
            Assert.assertArrayEquals("Numerical deviation above 1E-4", fidApodizedTRAFS, apodizedSpectrum.getFid(),1E-4);
            // consider removing this bit below
            boolean controlSTERR = false;
            int count = 0;
            for (int i =0; i< apodizedSpectrum.getFid().length;i++){
                data.add(i,apodizedSpectrum.getFid()[i]);
                if(Math.abs(apodizedSpectrum.getFid()[i]-fidApodizedTRAFS[i])>1E-5) {
                    if(!(controlSTERR) ){
                    System.err.println("TRAFS apodization with deviation larger than 1E-5");
                    System.err.println(i+" "+fidApodizedTRAFS[i]+" "+ apodizedSpectrum.getFid()[i]);
                    controlSTERR=true;
                    }
                    count++;
                }
            }
            System.out.println("TRAFS apodization deviations above 1E-5: "+count);
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        XYSeriesCollection dataset = new XYSeriesCollection(dataR);
        dataset.addSeries(data);
        JFreeChart chart = createChart(dataset);
        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart-trafs.png"), chart, 864 ,1152);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    private JFreeChart createChart(final XYDataset dataset) {
        final JFreeChart chart = ChartFactory.createScatterPlot(
                "Fid",                  // chart title
                "Hz?",                      // x axis label
                "Intensity",                      // y axis label
                dataset,                  // data
                PlotOrientation.HORIZONTAL.VERTICAL,
                true,                     // include legend
                true,                     // tooltips
                false                     // urls
        );

        XYPlot plot = (XYPlot) chart.getPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        renderer.setSeriesShapesVisible(0,false);


        plot.setRenderer(renderer);
        return chart;
    }

}
