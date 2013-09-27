/*
 * Copyright (c) 2013 EMBL, European Bioinformatics Institute.
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

package uk.ac.ebi.nmr.fid.tools;

import junit.framework.Assert;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;
import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Serializer;
import uk.ac.ebi.nmr.fid.Spectrum;

import java.io.File;
import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 24/05/2013
 * Time: 17:26
 * To change this template use File | Settings | File Templates.
 */
public class FastFourierTransform1DTest {


    /**
     * Test the fft on the putrescine experimental data, using the signal after the digital filter
     * @throws Exception
     */
    @Test
    public void testFFTBmse000109() throws Exception {


        double [] fidBmse000109DspCorrected= (double[]) Serializer.loadSerializedObject(
                FastFourierTransform1DTest.class.getClassLoader()
                .getResourceAsStream("data/bmse000109/1h/fid-DSPCorr.ser"));

        double [] fidBmse000109FFT1r= (double[]) Serializer.loadSerializedObject (
                FastFourierTransform1DTest.class.getClassLoader()
                .getResourceAsStream("data/bmse000109/1h/fid-fft-1r.ser"));
        double [] fidBmse000109FFT1i= (double[]) Serializer.loadSerializedObject (
                FastFourierTransform1DTest.class.getClassLoader()
                .getResourceAsStream("data/bmse000109/1h/fid-fft-1i.ser"));

        Acqu acquisition = new Acqu(Acqu.Spectrometer.BRUKER);
        acquisition.setAcquisitionMode(Acqu.AcquisitionMode.DISP);
        acquisition.setAquiredPoints(fidBmse000109DspCorrected.length);
        acquisition.setSpectralWidth(12.9911091032519);
        acquisition.setTransmiterFreq(499.842349248);
        acquisition.setDspGroupDelay(67.985595703125);// required for this DSP version
        acquisition.setDspFirmware(20);
        acquisition.setDspDecimation(3080);
        Spectrum spectrum = new Spectrum(fidBmse000109DspCorrected,acquisition);

        FastFourierTransform fft = new FastFourierTransform1D(spectrum);

        Spectrum transformedSpectrum = fft.computeFFT();


        XYSeries data = new XYSeries("FFT spectra");

        for (int i = 0 ; i< transformedSpectrum.getImaginaryChannelData().length ; i++)
            data.add(i,transformedSpectrum.getImaginaryChannelData()[i]);

        XYSeries dataR = new XYSeries("R data");

        double minDataR=0;
        double maxDataR=0;

        for (int i =0 ; i < fidBmse000109FFT1i.length ; i++){
            minDataR=Math.min(minDataR,fidBmse000109FFT1i[i]);
            minDataR=Math.min(minDataR, fidBmse000109FFT1r[i]);
            maxDataR=Math.max(maxDataR, fidBmse000109FFT1i[i]);
            maxDataR=Math.max(maxDataR, fidBmse000109FFT1r[i]);
        }

        for (int i = 0 ; i< fidBmse000109FFT1i.length ; i++) {
            dataR.add(i,fidBmse000109FFT1i[i]/(maxDataR-minDataR));
            fidBmse000109FFT1r[i]=fidBmse000109FFT1r[i]/(maxDataR-minDataR);
            fidBmse000109FFT1i[i]=fidBmse000109FFT1i[i]/(maxDataR-minDataR);
        }
        // calculate the rmsd between the current trasform and the one obtained from R
        Assert.assertTrue("Lengths of transfomed data differs",
                fidBmse000109FFT1r.length==transformedSpectrum.getRealChannelData().length);
        Assert.assertTrue("The transformed values differ above the limit",
                calculateRMSD(transformedSpectrum.getRealChannelData(),fidBmse000109FFT1r) < 0.02);
        Assert.assertTrue("Lengths of transfomed data differs",
                fidBmse000109FFT1i.length==transformedSpectrum.getImaginaryChannelData().length);
        Assert.assertTrue("The transformed values differ above the limit",
                calculateRMSD(transformedSpectrum.getImaginaryChannelData(),fidBmse000109FFT1i)< 0.02);

        XYSeriesCollection dataset = new XYSeriesCollection(data);
        dataset.addSeries(dataR);
        JFreeChart chart = createChart(dataset);
        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart-fft-ispectra.png"), chart, 864, 1152);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }

    /**
     * Test the fft on the putrescine experimental data, using the complete fid but setting the fft to use only
     * the signal after the digital filter.
     *
     * @throws Exception
     */
    @Test
    public void testFFTBmse000109Raw() throws Exception {


        double [] fidBmse000109DSPCorr= (double[]) Serializer.loadSerializedObject(
                FastFourierTransform1DTest.class.getClassLoader()
                .getResourceAsStream("data/bmse000109/1h/fid-DSPCorr-full.ser"));

        double [] fidBmse000109FFT1r= (double[]) Serializer.loadSerializedObject (
                FastFourierTransform1DTest.class.getClassLoader()
                .getResourceAsStream("data/bmse000109/1h/fid-fft-1r-full.ser"));
        double [] fidBmse000109FFT1i= (double[]) Serializer.loadSerializedObject (
                FastFourierTransform1DTest.class.getClassLoader()
                .getResourceAsStream("data/bmse000109/1h/fid-fft-1i-full.ser"));

        Acqu acquisition = new Acqu(Acqu.Spectrometer.BRUKER);
        acquisition.setAcquisitionMode(Acqu.AcquisitionMode.DISP);
        acquisition.setAquiredPoints(fidBmse000109DSPCorr.length);
        acquisition.setSpectralWidth(12.9911091032519);
        acquisition.setTransmiterFreq(499.842349248);
        acquisition.setDspGroupDelay(67.985595703125);// required for this DSP version
        acquisition.setDspFirmware(20);
        acquisition.setDspDecimation(3080);
        Spectrum spectrum = new Spectrum(fidBmse000109DSPCorr,acquisition);

        System.out.println(fidBmse000109DSPCorr.length);
        FastFourierTransform fft = new FastFourierTransform1D(spectrum);
        Spectrum transformedSpectrum = fft.computeFFT(140);


        XYSeries data = new XYSeries("FFT spectra");

        for (int i = 70 ; i< transformedSpectrum.getRealChannelData().length ; i++)
            data.add(i-70,transformedSpectrum.getRealChannelData()[i]);

        XYSeries dataR = new XYSeries("R data");

        double minDataR=0;
        double maxDataR=0;

        for (int i =0 ; i < fidBmse000109FFT1i.length ; i++){
            minDataR=Math.min(minDataR,fidBmse000109FFT1i[i]);
            minDataR=Math.min(minDataR, fidBmse000109FFT1r[i]);
            maxDataR=Math.max(maxDataR, fidBmse000109FFT1i[i]);
            maxDataR=Math.max(maxDataR, fidBmse000109FFT1r[i]);
        }

        for (int i = 0 ; i< fidBmse000109FFT1i.length ; i++) {
            dataR.add(i,fidBmse000109FFT1r[i]/(maxDataR-minDataR));
            fidBmse000109FFT1r[i]=fidBmse000109FFT1r[i]/(maxDataR-minDataR);
            fidBmse000109FFT1i[i]=fidBmse000109FFT1i[i]/(maxDataR-minDataR);
        }
        // calculate the rmsd between the current trasform and the one obtained from R
        Assert.assertTrue("Lengths of transfomed data differs",
                fidBmse000109FFT1r.length==(transformedSpectrum.getRealChannelData().length-70));
        Assert.assertTrue("The transformed values differ above the limit",
                calculateRMSD(transformedSpectrum.getRealChannelData(),fidBmse000109FFT1r) < 0.03);
        Assert.assertTrue("Lengths of transfomed data differs",
                fidBmse000109FFT1i.length==(transformedSpectrum.getImaginaryChannelData().length-70));
        Assert.assertTrue("The transformed values differ above the limit",
                calculateRMSD(transformedSpectrum.getImaginaryChannelData(),fidBmse000109FFT1i)< 0.03);
        System.out.println(fidBmse000109FFT1r.length+" "+transformedSpectrum.getRealChannelData().length);
        XYSeriesCollection dataset = new XYSeriesCollection(data);
        dataset.addSeries(dataR);
        JFreeChart chart = createChart(dataset);

        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart-fft-spectra.png"), chart, 864, 1152);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
    private double calculateRMSD(double[] realChannelData, double[] fidFFTr) {
        double rmsd=0;
        // if only part of the fid is transformed, case of offset
        int shortestLength=Math.min(fidFFTr.length,realChannelData.length);
        for (int i =0 ; i < shortestLength; i++)
            rmsd+= Math.pow(fidFFTr[i] - realChannelData[i+(realChannelData.length-shortestLength)],2);

        return Math.sqrt(rmsd/shortestLength);
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
