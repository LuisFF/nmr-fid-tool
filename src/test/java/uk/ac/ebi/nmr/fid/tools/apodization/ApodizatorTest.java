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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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

    static double dw=0.01;

    XYSeries data;
    static XYSeries dataR;

    @BeforeClass
    public static void loadExternalData () {

        final Pattern DATAROW_PATTERN = Pattern.compile("(-?\\d+\\.?\\d*e?-?\\+?\\d*)");
        final Pattern DATA_PATTERN = Pattern
                .compile("(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)"+
                "\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)"+
                "\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)");
        
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(ApodizatorTest.class.getClassLoader()
                .getResourceAsStream("data/fid-simulated.csv")));
        int lines = 0;
        // count the number of lines
        // which corresponds to the number of data points if one removes the header line
        try {

            bufferedReader.readLine();
            while (bufferedReader.readLine() != null) lines++;

            bufferedReader.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        fid = new double[lines];
        double[] fid4spectrum = new double[lines];
        fidApodizedEM = new double[lines];
        fidApodizedGM = new double[lines];
        fidApodizedLGM = new double[lines];
        fidApodizedTRAF = new double[lines];
        fidApodizedTRAFS = new double[lines];
        realFTObserved = new double[lines];
        imaginaryFTObserved = new double[lines];

        bufferedReader = new BufferedReader(new InputStreamReader(ApodizatorTest.class.getClassLoader()
                .getResourceAsStream("data/fid-simulated.csv")));



        try {
            String line = bufferedReader.readLine();
            int index = 0;
            Matcher matcher = null;
            dataR = new XYSeries("dataR");

            while (line != null){
                matcher = DATAROW_PATTERN.matcher(line);
                if(matcher.find()){
                    int matchEnd=0;
                    int times =0;
                    while (matcher.find(matchEnd)){
                        matchEnd=matcher.end();
//                        System.out.println("Pattern found "+ new String(new char[times]).replace("\0","\t" )
//                                +Double.parseDouble(matcher.group(1)));
                        switch (times){
                            case 0 :
                                fid[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 1 :
                                fid4spectrum[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 2 :
                                realFTObserved[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 3 :
                                imaginaryFTObserved[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 4 :
                                fidApodizedEM[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 5 :
                                fidApodizedGM[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 6 :
                                fidApodizedLGM[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 7 :
                                fidApodizedTRAF[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 8 :
                                fidApodizedTRAFS[index]=Double.parseDouble(matcher.group(1));
                                break;
                            default:
                                break;
                        }
                        ++times;
                    }
                    dataR.add(index*dw,fid4spectrum[index]);
                    ++index;
                }

                line=bufferedReader.readLine();
            }
            System.out.println(index);
            acquisition = new Acqu(Acqu.Spectrometer.BRUKER);
            // The test data I have is sequencial, meaning that I have only the FID of the real channel
            // The simultaneous would have the value of the real channel followed by the corresponding value of the
            // imaginary chanel
            acquisition.setAcquisitionMode(Acqu.AcquisitionMode.SEQUENTIAL);
            System.out.println(acquisition.getAcquisitionMode().toString());
            acquisition.setAquiredPoints(fid.length);
            acquisition.setTransmiterFreq(8 + 1 / 3); // for bmse000109 is app. 500
            acquisition.setSpectralWidth(6.25); // for bmse000109 is app. 13
            processing = new Proc(acquisition);
            processing.setLineBroadening(1);
            processing.setGbFactor(0.03);
            System.out.println("DW: " +processing.getDwellTime());
            System.out.println("LB: " +processing.getLineBroadening());
            spectrum = new Spectrum(fid4spectrum,acquisition);

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        }

        @Test
        public void TestExponentialApodizatorSequential(){
            data = new XYSeries("data");

            Apodizator expApodizator= new ExponentialApodizator(spectrum);


            try {
                Spectrum apodizedSpectrum = expApodizator.calculate(1);
            boolean controlSTERR = false;
            for (int i =0; i< apodizedSpectrum.getFid().length;i++){
                data.add(i*dw,apodizedSpectrum.getFid()[i]);
                if(!(controlSTERR) && Math.abs(apodizedSpectrum.getFid()[i]-fidApodizedEM[i])>1E-12) {
                    System.err.println("Exponential apodization with deviation larger than 1E-12");
                    System.err.println(i+" "+ fidApodizedEM[i]+" "+ apodizedSpectrum.getFid()[i]);
                    controlSTERR=true;
                }
                Assert.assertTrue("Numerical deviation above 1E-12",
                        Math.abs(apodizedSpectrum.getFid()[i]-fidApodizedEM[i])<1E-12);

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
    public void TestGaussianApodizatorSequential(){
        data = new XYSeries("data");
        Apodizator gaussApodizator= new GaussianApodizator(spectrum);

        try {
            Spectrum apodizedSpectrum = gaussApodizator.calculate(1);
            boolean controlSTERR = false;
            for (int i =0; i< apodizedSpectrum.getFid().length;i++){
                data.add(i*dw,apodizedSpectrum.getFid()[i]);
                if(Math.abs(apodizedSpectrum.getFid()[i]-fidApodizedGM[i])>1E-12) {
                    if(!controlSTERR){
                        System.err.println("Gaussian apodization with deviation larger than 1E-12");
                        System.err.println(i+" "+ fidApodizedGM[i]+" "+ apodizedSpectrum.getFid()[i]);
                        controlSTERR=true;
                    }
                }
                Assert.assertTrue("Numerical deviation above 1E-12",
                        Math.abs(apodizedSpectrum.getFid()[i] - fidApodizedGM[i]) < 1E-12);

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
    public void TestLorentzGausApodizatorSequential(){
        data = new XYSeries("data");
        Apodizator lgApodizator= new LorentzGausApodizatior(spectrum);

        try {
            Spectrum apodizedSpectrum = lgApodizator.calculate(-0.3);
            boolean controlSTERR = false;
            int count=0;
            for (int i =0; i< apodizedSpectrum.getFid().length;i++){
                data.add(i*dw,apodizedSpectrum.getFid()[i]);
                if(Math.abs(apodizedSpectrum.getFid()[i]-fidApodizedLGM[i])>1E-3) {
                    if(!(controlSTERR)){
                    System.err.println("Lorentz-Gauss apodization with deviation larger than 1E-3");
                    System.err.println(i+" "+fidApodizedLGM[i]+" "+ apodizedSpectrum.getFid()[i]);
                        controlSTERR=true;
                    }
                    count++;
                }
            }
            System.out.println("Lorentz-Gauss apodization deviations above 1E-3: "+count);
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
            boolean controlSTERR = false;
            int count = 0;
            for (int i =0; i< apodizedSpectrum.getFid().length;i++){
                data.add(i*dw,apodizedSpectrum.getFid()[i]);
                Assert.assertTrue("Numerical deviation above 1E-3",
                        Math.abs(apodizedSpectrum.getFid()[i] - fidApodizedTRAF[i]) < 1E-3);
                if(Math.abs(apodizedSpectrum.getFid()[i]-fidApodizedTRAF[i])>1E-4) {
                    if(!(controlSTERR)){
                    System.err.println("TRAF apodization with deviation larger than 1E-4");
                    System.err.println(i+" "+fidApodizedTRAF[i]+" "+ apodizedSpectrum.getFid()[i]);
                    controlSTERR=true;
                    }
                    count++;
                }
            }
            System.out.println("TRAF apodization deviations above 1E-4: "+count);
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
    public void TestTRAFSApodizatorSequential(){
        data = new XYSeries("data");
        Apodizator tarfsApodizator= new TrafsApodizator(spectrum);

        try {
            Spectrum apodizedSpectrum = tarfsApodizator.calculate(0.06);
            boolean controlSTERR = false;
            int count = 0;
            for (int i =0; i< apodizedSpectrum.getFid().length;i++){
                data.add(i*dw,apodizedSpectrum.getFid()[i]);
                Assert.assertTrue("Numerical deviation above 1E-3",
                        Math.abs(apodizedSpectrum.getFid()[i] - fidApodizedTRAFS[i ]) < 1E-3);
                if(Math.abs(apodizedSpectrum.getFid()[i]-fidApodizedTRAFS[i])>1E-4) {
                    if(!(controlSTERR) ){
                    System.err.println("TRAFS apodization with deviation larger than 1E-4");
                    System.err.println(i+" "+fidApodizedTRAFS[i]+" "+ apodizedSpectrum.getFid()[i]);
                    controlSTERR=true;
                    }
                    count++;
                }
            }
            System.out.println("TRAFS apodization deviations above 1E-4: "+count);
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
