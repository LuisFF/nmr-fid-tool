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


import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Max;
import org.apache.commons.math3.stat.descriptive.rank.Min;
import uk.ac.ebi.nmr.fid.Spectrum;

import java.util.ArrayList;
import java.util.List;

/**
 * Implementation of a baseline corrector based on the method from Golotvin and Williams (2000)
 *
 * @author  Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 29/05/2013
 * Time: 15:48
 * To change this template use File | Settings | File Templates.
 */
public class GolotvinWilliamsBaselineCorrector implements BaselineCorrector {


    @Override
    public Spectrum correctBaseline(Spectrum spectrum) {

        boolean [] baselineIndexes = baslineDetection(spectrum.getRealChannelData());

        double [] baselineModel = calculateBaselineModel(baselineIndexes,spectrum.getRealChannelData());

        double [] subtractedReal = new double [spectrum.getRealChannelData().length];
        System.arraycopy(spectrum.getRealChannelData(),0,subtractedReal,0,spectrum.getRealChannelData().length);

        for (int i = 0 ; i < spectrum.getRealChannelData().length;i++){
            subtractedReal[i]-=baselineModel[i];
        }

        spectrum.setBaseline(baselineIndexes);
        spectrum.setBaselineModel(baselineModel);
        spectrum.setRealChannelData(subtractedReal);
        return spectrum;  //To change body of implemented methods use File | Settings | File Templates.
    }

    private double[] calculateBaselineModel(boolean[] baselineIndexes, double[] realChannelData) {
        int windowSize = 41;
        double [] convolutionWindow = new double[windowSize];
        double [] baselineModel= new double[realChannelData.length];
        LinearInterpolator interpolator = new LinearInterpolator();
        double [] indexes= new double[2];
        indexes[0]=-1;
        double [] intesities= new double[2];
        List<Double> window = new ArrayList<Double>(windowSize);
        double smoothValue = 0;
        int lastVisitedIndex = 0;

        // preload the window list with points at the end of the spectrum
        for(int i = realChannelData.length-1; i> 0 ; i--){
            if(baselineIndexes[i]){
                // is this going to work????
                window.add(0,realChannelData[i]/windowSize);
                smoothValue+=realChannelData[i]/windowSize;
            }
            if(window.size()>=(windowSize/2-1))
                break;
        }

        for (int i = 0; i<realChannelData.length; i++){
            // for each baseline point calculate the smoothed spectrum
            if(baselineIndexes[i]){
                // first baseline point considered so the window is not fully loaded
                if(window.size()<windowSize){
                    window.add(realChannelData[i]/windowSize);
                    smoothValue+=realChannelData[i]/windowSize;
                    // load the window with the right side values
                    for(int j = i+1; j<realChannelData.length; j++){
                        if(baselineIndexes[j]){
                            lastVisitedIndex=j;
                            window.add(realChannelData[j]/windowSize);
                            smoothValue+=realChannelData[j]/windowSize;
                        }
                        if(window.size()>=windowSize)
                            break;
                    }
                    if(lastVisitedIndex==realChannelData.length)
                        System.err.println("Problem loading the window in the baseline correction");
                } else {
                    // just update the window and the smoothvalue
                    smoothValue-=window.get(0);
                    window.remove(0);
                    while(lastVisitedIndex < (realChannelData.length-1)){
                        lastVisitedIndex++;
                        if(baselineIndexes[lastVisitedIndex]){
                            window.add(realChannelData[lastVisitedIndex]/windowSize);
                            smoothValue+=realChannelData[lastVisitedIndex]/windowSize;
                        }
                        if(window.size()==windowSize)
                            break;
                    }
                    // if we were not able to fill in the window with values
                    if(window.size()<windowSize){
                        // this should only happen if we reach the edge of the spectrum...
                        lastVisitedIndex=-1;
                        while(lastVisitedIndex < realChannelData.length){
                            lastVisitedIndex++;
                            if(baselineIndexes[lastVisitedIndex]){
                                window.add(realChannelData[lastVisitedIndex]/windowSize);
                                smoothValue+=realChannelData[lastVisitedIndex]/windowSize;
                            }
                            if(window.size()==windowSize)
                                break;
                        }
                    }
                }
                baselineModel[i]=smoothValue;
                // check if we are coming from a non-baseline point and do a linear interpolation
                // for the non-baseline points in between
                // TODO !!!make sure this is not the first point!!!! but the first point would break this...
                if(i>0){
                if(!baselineIndexes[i-1] && indexes[0]>0){
                    indexes[1]=i;
                    intesities[1]=baselineModel[i];
                    PolynomialSplineFunction splineFunction = interpolator.interpolate(indexes,intesities);
                    for(int j = (int) indexes[0]; j< indexes[1] ;j++)
                        baselineModel[j]=splineFunction.value(j);
                }
                }
                indexes[0]=i;
                intesities[0]=baselineModel[i];
            }

        }
        // I still need to take care of the cases where edges of the spectrum are not baseline points...
        if(!baselineIndexes[0] || !baselineIndexes[baselineIndexes.length-1]){
            // fetch the of the extreme baseline indexes
            for(int i = baselineModel.length-1; i>0; i--){
                if(baselineIndexes[i]){
                    indexes[0]=i-baselineModel.length;
                    intesities[0]=baselineModel[i];
                    break;
                }
            }
            for(int i = 0 ; i< baselineModel.length;i++){
                if(baselineIndexes[i]){
                    indexes[1]=i;
                    intesities[1]=baselineModel[i];
                    break;
                }
            }
            PolynomialSplineFunction splineFunction = interpolator.interpolate(indexes,intesities);
            for(int j = (int) indexes[0]+1; j < 0 ;j++)
                baselineModel[j+baselineModel.length]=splineFunction.value(j);
            for(int j = 0; j <indexes[1];j++){
                baselineModel[j]=splineFunction.value(j);
            }

        }

        return baselineModel;  //To change body of created methods use File | Settings | File Templates.
    }

    private boolean[] baslineDetection(double[] realChannelData) {
        int rectangleWidth=60;
        int factor = 4;
        double noiseStandardDeviation = calculateSpectralNoise(realChannelData);
        boolean [] baseline = new boolean[realChannelData.length];
        Max maximum = new Max();
        Min minimum = new Min();
        for(int i = 0; i<realChannelData.length; i++){
            double [] datapoints = new double[rectangleWidth+1];
            //Let's assume that the spectrum is a continuum of points like in Friedrichs, M. (1995)
            if(i<(rectangleWidth/2)){
                System.arraycopy(realChannelData,0,datapoints,0,i+rectangleWidth/2+1);
                // take the missing points from the end of the spectra
                System.arraycopy(realChannelData,realChannelData.length-(rectangleWidth/2-i),
                        datapoints,i+rectangleWidth/2+1,rectangleWidth-(i+rectangleWidth/2));
            } else if((realChannelData.length-i)<(rectangleWidth/2+1)){
                //collect the last i points plus the the left side points of i
                System.arraycopy(realChannelData,i-rectangleWidth/2,datapoints,0,realChannelData.length-i+rectangleWidth/2);
                // take the missing points from the beginning of the spectra
                System.arraycopy(realChannelData,0,datapoints,
                        realChannelData.length-i+rectangleWidth/2,rectangleWidth/2-(realChannelData.length-i)+1);
            } else {
                System.arraycopy(realChannelData,i-rectangleWidth/2,datapoints,0,rectangleWidth+1);
            }
            baseline[i]=(maximum.evaluate(datapoints)-minimum.evaluate(datapoints))<(factor*noiseStandardDeviation);
        }

        return baseline;  //To change body of created methods use File | Settings | File Templates.
    }

    private double calculateSpectralNoise(double[] realChannelData) {
        int numberOfRegions=1;
        // define the window size starting with 32 datapoints
        // and shrinking it until the exact division is achieved
        for (int i =32; i > 0 ; i--){
            // check if the datapoints can be exactly divided into i regions
            if(realChannelData.length % i ==0 ){
                numberOfRegions = i;
                break;
            }
        }
        //set the noise to the maximum intensity
        Max max = new Max();
        StandardDeviation standardDeviation = new StandardDeviation();

        double noiseStandardDeviation=max.evaluate(realChannelData);

        //determine the standard deviation for each window and record the lowest standard deviation
        for (int i = 0; i<realChannelData.length; i+=realChannelData.length/numberOfRegions){
            double [] region = new double[realChannelData.length/numberOfRegions];
            System.arraycopy(realChannelData,i,region,0,realChannelData.length/numberOfRegions);
            noiseStandardDeviation=(noiseStandardDeviation<standardDeviation.evaluate(region))?
                    noiseStandardDeviation:standardDeviation.evaluate(region);
        }
        return noiseStandardDeviation;
    }


}
