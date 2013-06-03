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

package uk.ac.ebi.nmr.fid.tools.baseline;

import org.apache.commons.math3.stat.descriptive.rank.Max;
import org.apache.commons.math3.stat.descriptive.rank.Min;
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
import org.junit.Test;
import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Serializer;
import uk.ac.ebi.nmr.fid.Spectrum;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * Test classes for the baseline corrector
 *
 * @author  Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 03/06/2013
 * Time: 12:26
 * To change this template use File | Settings | File Templates.
 */
public class GolotvinWilliamsBaselineCorrectorTest {

    @Test
    public void testCorrectBaseline() throws Exception {

        double [] fidBmse000109DSPCorr= (double[]) Serializer.loadSerializedObject(
                GolotvinWilliamsBaselineCorrectorTest.class.getClassLoader()
                        .getResourceAsStream("data/bmse000109/1h/fid-DSPCorr-full.ser"));

        double [] fidBmse000109FFT1r= (double[]) Serializer.loadSerializedObject(
                GolotvinWilliamsBaselineCorrectorTest.class.getClassLoader()
                        .getResourceAsStream("data/bmse000109/1h/fid-fft-phased-1r.ser"));

        // I used the normalized spectra in R
        double maximum = new Max().evaluate(fidBmse000109FFT1r);
        double minimum = new Min().evaluate(fidBmse000109FFT1r);
        for(int i=0; i < fidBmse000109FFT1r.length;i++)
            fidBmse000109FFT1r[i]=fidBmse000109FFT1r[i]/(maximum-minimum);

        double [] fidBmse000109FFT1i= (double[]) Serializer.loadSerializedObject(
                GolotvinWilliamsBaselineCorrectorTest.class.getClassLoader()
                        .getResourceAsStream("data/bmse000109/1h/fid-fft-phased-1i.ser"));

        double [] fidBmse000109FFT1rBaselineCorrected= (double[]) Serializer.loadSerializedObject(
                GolotvinWilliamsBaselineCorrectorTest.class.getClassLoader()
                        .getResourceAsStream("data/bmse000109/1h/fid-fft-1r-baselinecorrected.ser"));

        double [] fidBmse000109FFT1rBaselineModel= (double[]) Serializer.loadSerializedObject(
                GolotvinWilliamsBaselineCorrectorTest.class.getClassLoader()
                        .getResourceAsStream("data/bmse000109/1h/fid-fft-1r-baselinemodel.ser"));

        boolean [] fidBmse000109FFT1rBaselineIndexes= (boolean[]) Serializer.loadSerializedObject(
                GolotvinWilliamsBaselineCorrectorTest.class.getClassLoader()
                        .getResourceAsStream("data/bmse000109/1h/fid-fft-1r-baseline.ser"));

        Acqu acquisition = new Acqu(Acqu.Spectrometer.BRUKER);
        acquisition.setAcquisitionMode(Acqu.AcquisitionMode.DISP);
        acquisition.setAquiredPoints(fidBmse000109DSPCorr.length);
        acquisition.setSpectralWidth(12.9911091032519);
        acquisition.setTransmiterFreq(499.842349248);
        acquisition.setDspGroupDelay(67.985595703125);// required for this DSP version
        acquisition.setDspFirmware(20);
        acquisition.setDspDecimation(3080);
        Spectrum spectrum = new Spectrum(fidBmse000109DSPCorr,acquisition);
        spectrum.setRealChannelData(fidBmse000109FFT1r);
        spectrum.setImaginaryChannelData(fidBmse000109FFT1i);
        
        BaselineCorrector baselineCorrector = new GolotvinWilliamsBaselineCorrector();
        
        Spectrum baselineCorrectedSpectrum = baselineCorrector.correctBaseline(spectrum);

        Assert.assertTrue("Baseline points differently defined",Arrays.equals(fidBmse000109FFT1rBaselineIndexes,
                baselineCorrectedSpectrum.getBaseline()));

        Assert.assertArrayEquals("Baseline model has a deviation larger than 1E-4", fidBmse000109FFT1rBaselineModel,
                baselineCorrectedSpectrum.getBaselineModel(),1E-4);

        Assert.assertArrayEquals("Corrected data has a deviation larger than 1E-4", fidBmse000109FFT1rBaselineCorrected,
                baselineCorrectedSpectrum.getRealChannelData(),1E-4);

        XYSeries data = new XYSeries("baseline corrected");
        XYSeries dataR = new XYSeries("R data");

        for (int i = 0 ; i<baselineCorrectedSpectrum.getRealChannelData().length; i++){
            data.add(i,baselineCorrectedSpectrum.getRealChannelData()[i]);
        }
        for(int i = 0 ; i< fidBmse000109FFT1r.length;i++){
            dataR.add(i,fidBmse000109FFT1r[i]);
        }

        XYSeriesCollection dataset = new XYSeriesCollection(data);
        dataset.addSeries(dataR);
        JFreeChart chart = createChart(dataset);

        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart-baseline-spectra.png"), chart, 864, 1152);
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
