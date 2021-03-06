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

package uk.ac.ebi.nmr.execs;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import edu.uchc.connjur.core.AnalysisDomain;
import edu.uchc.connjur.core.DimInfo;
import edu.uchc.connjur.spectrumtranslator.DataControl;
import edu.uchc.connjur.spectrumtranslator.SwDataSet;
import edu.uchc.connjur.spectrumtranslator.bruker.BrukerDataSetReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;
import uk.ac.ebi.nmr.fid.Spectrum;
import uk.ac.ebi.nmr.fid.io.AcquReader;
import uk.ac.ebi.nmr.fid.io.BrukerAcquReader;
import uk.ac.ebi.nmr.fid.io.ConnjurFidReader;
import uk.ac.ebi.nmr.fid.io.FidReader;

import javax.swing.*;
import java.io.File;
import java.io.InputStream;
import java.util.Iterator;

/**
 * scratch class to vizualize and do some experiments with the data
 *
 * @author Luis F. de Figueiredo
 */

public class PlotFID02 extends JFrame {

    private XYSeriesCollection dataset;

    public static void main(String[] args) throws Exception {
        new PlotFID02().showGraph();
    }

    public PlotFID02() throws Exception {

        dataset = new XYSeriesCollection();



//        initialApproach();
        connjurIntegration();


    }

    private void connjurIntegration() throws Exception {
        final XYSeries data = new XYSeries("data");
        final XYSeries data2 = new XYSeries("pivot");
        BrukerDataSetReader dataSetReader = new BrukerDataSetReader();

        dataSetReader.setDataSource(new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/" +
                "resources/examples/file_formats/bmse000109/1H/"));
//        dataSetReader.setDataSource(new File("/Users/ldpf/Data/Databases/bmrb/tar_01/4_aminobenzoic_acid/nmr/bmse000066/1H/"));
//        dataSetReader.setDataSource(new File("/Users/ldpf/Data/EIPOD/nmr/purecompound-201204-05/nmr/lf_gliotoxin_600/1/"));


        DataControl dataControl = new DataControl();
        SwDataSet dataSet = dataSetReader.read(null,dataControl);
//        SwDataSet dataSet = dataSetReader.read(true);

        System.out.println("Number of axes " + dataSet.getNumberAxes());

        System.out.println("Data type if the binary file: " + dataSet.getSourceDataType());
        DimInfo dimInfo = dataSet.getDataByDim(0);
        System.out.println("First order Phase :" + dimInfo.getFirstOrderPhase());
        System.out.println("Zero order Phase :" + dimInfo.getZeroOrderPhase());
        System.out.println("Total number of values (Re + Im):" + dimInfo.getTotalPoints());


        // get the frequency domain
        boolean isFrequency = false;
        for (int i = 0; i < dataSet.getNumberAxes() && !isFrequency; i++) {
            if (dataSet.getDataByDim(i).getAnalysisDomain() == AnalysisDomain.FREQUENCY) {
                isFrequency = true;
            }
        }
        // the real channel is accessible in the first fid?
        // the imaginary channel is accessible in the second fid?
        Iterator<edu.uchc.connjur.core.Fid> fidIterator = dataSet.fidIterator();
        if (fidIterator.hasNext()){
            // get the First fid - there should be more Fid for 2 dimensional spectra

            edu.uchc.connjur.core.Fid fidReal = fidIterator.next();
            //prepare for the appodization

            Acqu acquisition = new Acqu(Acqu.Spectrometer.BRUKER);
            acquisition.setAcquisitionMode(Acqu.AcquisitionMode.SIMULTANIOUS);
            System.out.println(acquisition.getAcquisitionMode().toString());
            acquisition.setAquiredPoints((fidReal.getValues().length-140)*2);
            acquisition.setTransmiterFreq(500); // for bmse000109 is app. 500
            acquisition.setSpectralWidth(13); // for bmse000109 is app. 13
            Proc processing = new Proc(acquisition);
            processing.setLineBroadening(1);
            System.out.println("DW: " +processing.getDwellTime());
            System.out.println("LB: " +processing.getLineBroadening());


            int offSet = 141; // number of points ignored in the left side of the fid
            System.out.println("number of points in the re channel: " + (fidReal.getValues().length));
            edu.uchc.connjur.core.Fid fidIm = fidIterator.next();
            System.out.println("number of points in the im channel: "+ (fidIm.getValues().length));
            double  [] fidTransformed= new double[(fidReal.getSize()-140)*2];
//            float [] fidImValues= new float[fidIm.getSize()-offSet];
            // interlace the real and imaginary points into a single array
            for (int i = 140; i < fidReal.getSize(); i++ ){
                fidTransformed[(i-140)*2]=fidReal.getValue(i);// insert the real (odd entries)
                fidTransformed[(i-140)*2+1]=fidIm.getValue(i);// insert the imaginary (even entries)
            }

//            ExponentialApodizator apodizator = new ExponentialApodizator(fidTransformed,processing);


            DoubleFFT_1D fftd = new DoubleFFT_1D(fidReal.getValues().length-140);
//            double[] spectrumRaw = apodizator.calculate(1.0);
            double[] spectrumRaw = fidTransformed;
            fftd.complexInverse(spectrumRaw, true);
            // reorganize the spectrum
            // swapt the left with the rigth side and invert the right side
            // do this at the global scale so that it can be used for the phase correction
            //TODO find the exact place to split the spectrum
//fftMirrored[1:((length(fftValue)+1)/2)]=fftValue[((length(fftValue)+1)/2):length(fftValue)];
//fftMirrored[((length(fftValue)+1)/2):length(fftValue)]=fftValue[1:((length(fftValue)+1)/2)];

            AcquReader acquReader = new BrukerAcquReader(new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/" +
                    "resources/examples/file_formats/bmse000109/1H/acqus"));
            Acqu realAcqu = acquReader.read();
            System.out.println(realAcqu.getSpectralWidth());
            System.out.println("Total point of fid transformed (re+im)"+fidTransformed.length);
            // get only the real values, which should be in the odd entries

            // hmmm
            double teta0 = -0.33*2*Math.PI;
            double teta1 = 2.3*Math.PI;
            int pivot = 3850;
            double [] spectrum = phasecorrection(spectrumRaw,teta0, teta1, pivot);
            double[] spectrumMirrowed= new double[spectrum.length];
            System.out.println(spectrum.length/2);
            if(pivot<spectrum.length/2)
                data2.add((float)(spectrum.length/2-pivot)/(spectrum.length)*realAcqu.getSpectralWidth()-1.3,-3000);
            else
                data2.add((float)((-spectrum.length/2 + pivot)/(spectrum.length))*realAcqu.getSpectralWidth()-1.3,-3000);
            // swap sides
            for (int i = 0; i< spectrum.length/2; i++){ // check how to do with odd lengths...
                spectrumMirrowed[spectrum.length/2+i]=spectrum[i];
                spectrumMirrowed[i]=spectrum[spectrum.length/2+i];
            }


//            fftd.realInverse(fidImValues,true);
            // create the quadrature values for the spectra
//            quadLeft=Re(realPart)-Im(imagPart);
//            quadRight=Re(realPart)+Im(imagPart);






            // do the appodization
//            ApodizationTool apodizationTool = new ApodizationTool(convertFloatsToDoubles(fidRealValues),
//                    Acqu.AcquisitionMode.SIMULTANIOUS,processing);
//            double [] appodizedValues =apodizationTool.exponential();
//            double [] appodizedValues =convertFloatsToDoubles(fidRealValues);
            // prepare for the FFT
//            DoubleFFT_1D dfftd = new DoubleFFT_1D(appodizedValues.length/2);
//            dfftd.complexForward(appodizedValues);

            // do the phase correction
//            PhaseCorrectionTool phaseCorrectionTool = new PhaseCorrectionTool(appodizedValues, new Proc(acquisition));
//            double [] spectrumPhased = phaseCorrectionTool.zeroOrderPhasing(dimInfo.getZeroOrderPhase());
//            double [] spectrumPhased = phaseCorrectionTool.firstOrderPhasing(dimInfo.getFirstOrderPhase());
            // get just the real part of the points
//            double[] realPart = new double[appodizedValues.length/2];
//            int realIndex = 0;
//            for (int i = 0; i < appodizedValues.length-1; i+=2) {
//                realPart[realIndex++] = appodizedValues[i];
//            }
//            // put everything in data
//            for (int i=realPart.length-1; i>=0 ; i--){
//                data.add(realPart.length-i,realPart[i]);
//            }


//            putrescine case
//            double teta0 = -0.08*2*Math.PI;
//            double teta1 = -2.0*Math.PI;
//            int pivot = 6400;
            //            putrescine case - Jacobsen Approach
            // 2 peaks good
//            double teta0 = 0.05*2*Math.PI;
//            double teta1 = -0.5*Math.PI;
//            int pivot = 8650;
            // 2other peaks phased
//            double teta0 = 0.05*2*Math.PI;
//            double teta1 = -2.4*Math.PI;
//            int pivot = 8650;
//            // hmmm
//            double teta0 = -0.45*2*Math.PI;
//            double teta1 = 0*Math.PI;
//            int pivot = 8650;
            //4_aminobenzoic_acid
//            double teta0 = 0.1*2*Math.PI;
//            double teta1 = -2.5*Math.PI;
//            int pivot = 1500;
//            double [] spectrum = phasecorrection(spectrumMirrowed,teta0, teta1, pivot);
//            data2.add((float)(spectrum.length - pivot)/(spectrum.length)*realAcqu.getSpectralWidth()-1.3,-3000);

            //TODO break the spectra into half and mirror the left side...
            // this works more or less ok
//            for (int i = 0; i< spectrumRaw.length/2; i++){
//                data.add((float)(spectrumRaw.length/2 - i)/(spectrumRaw.length/2)*realAcqu.getSpectralWidth()-1.65,spectrumRaw[i*2+1]);
////                data.add(i,spectrumRaw[i*2]);
//            }
            // second approach
            for (int i = 0; i< spectrumMirrowed.length; i++){
                //putrescine
                data.add((float)(spectrumMirrowed.length - i)/(spectrumMirrowed.length)*realAcqu.getSpectralWidth()-1.65,spectrumMirrowed[i]);
                // 4_aminobenzoic_acid
//                data.add((float)(spectrum.length - i)/(spectrum.length)*realAcqu.getSpectralWidth()-1.3,-spectrum[i]);
            }
        }


//            data.add(3, 2); //Point 1
//            data.add(1, 1); //Point 2
//            data.add(4, 1); //Point 3
//            data.add(2, 2); //Point 4

        dataset.addSeries(data);
        dataset.addSeries(data2);
    }

    private double[] phasecorrection(double[] spectrumRaw, double teta0, double teta1, int pivot) {
        double []  spectrum = new double [spectrumRaw.length/2];
        System.out.println(teta0+teta1*(0-pivot)/spectrum.length);
        for(int i =0 ; i < spectrum.length; i++){
            spectrum[i]=spectrumRaw[i*2]*Math.cos(teta0+teta1*(i-pivot)/spectrum.length)+
                    spectrumRaw[i*2+1]*Math.sin(teta0 + teta1 * (i - pivot) / spectrum.length);
        }
        return spectrum;
    }

    private void initialApproach() throws Exception {

        final XYSeries data = new XYSeries("data");
        Acqu acquisition = (new BrukerAcquReader(PlotFID02.class.getResourceAsStream("/examples/file_formats/bmse000109/1H/acqu"))).read();
        InputStream fidInput = PlotFID02.class.getResourceAsStream("/examples/file_formats/bmse000109/1H/fid");

        FidReader fidReader = new ConnjurFidReader(new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/"+
                "resources/examples/file_formats/bmse000109/1H/"), acquisition);
//        Fid fid = new FidReader(fidInput, acquisition).read();
        Spectrum spectrum = fidReader.read();
        System.out.println("numeber of points: " + spectrum.getFid().length);


        dataset.addSeries(data);
    }

    public void showGraph() {



        final JFreeChart chart = createChart(dataset);
        final ChartPanel chartPanel = new ChartPanel(chart);
        // invert axis
        chartPanel.getChart().getXYPlot().getDomainAxis().setInverted(true);
        // define the axis
//        chartPanel.getChart().getXYPlot().getRangeAxis().setRange(-10,50);
        chartPanel.getChart().getXYPlot().getRangeAxis().setVisible(false);
        chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
        final ApplicationFrame frame = new ApplicationFrame("Title");

        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    public JFreeChart createChart(final XYDataset dataset) {
        final JFreeChart chart = ChartFactory.createScatterPlot(
                "putrescine 1H spectra", // chart title
                "ppm", // x axis label
                "", // y axis label
                dataset, // data
                PlotOrientation.HORIZONTAL.VERTICAL,
                true, // include legend
                true, // tooltips
                false // urls
                );

        XYPlot plot = (XYPlot) chart.getPlot();

        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        renderer.setSeriesShapesVisible(0, false);


        plot.setRenderer(renderer);
        return chart;
    }
    public static double[] convertFloatsToDoubles(float[] input)
    {
        if (input == null)
        {
            return null; // Or throw an exception - your choice
        }
        double[] output = new double[input.length];
        for (int i = 0; i < input.length; i++)
        {
            output[i] = input[i];
        }
        return output;
    }
}
