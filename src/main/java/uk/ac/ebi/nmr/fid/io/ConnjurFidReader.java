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

package uk.ac.ebi.nmr.fid.io;

import edu.uchc.connjur.spectrumtranslator.DataControl;
import edu.uchc.connjur.spectrumtranslator.SwDataSet;
import edu.uchc.connjur.spectrumtranslator.bruker.BrukerDataSetReader;
import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Spectrum;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/**
 * Connjur integration to read the fid data
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 21/03/2013
 * Time: 12:00
 *
 */
public class ConnjurFidReader implements FidReader {


    private BrukerDataSetReader dataSetReader;
    private Acqu acquisition;

    public ConnjurFidReader(File file, Acqu acqu) {
        this.dataSetReader = new BrukerDataSetReader();
        this.acquisition=acqu;
        try {
            dataSetReader.setDataSource(file);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }

    @Override
    public Spectrum read() throws Exception {
        


        DataControl dataControl = new DataControl();
        SwDataSet dataSet = dataSetReader.read(null,dataControl);
        // extract the real and imaginary values of the FID
        Iterator<edu.uchc.connjur.core.Fid> fidIterator = dataSet.fidIterator();
        edu.uchc.connjur.core.Fid fidReal = fidIterator.next();
        edu.uchc.connjur.core.Fid fidIm = fidIterator.next();

        //put them together because the FID only accepts a list of alternating real + imaginary numbers
        double [] fidCombined = new double[fidReal.getSize()+fidIm.getSize()];
        for(int i = 0; i < fidReal.getSize();i++){
            fidCombined[i*2]= fidReal.getValue(i);
            fidCombined[i*2+1]= fidIm.getValue(i);
        }

        return new Spectrum(fidCombined, acquisition);
    }
}
