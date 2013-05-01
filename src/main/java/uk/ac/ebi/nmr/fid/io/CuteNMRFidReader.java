/**
 * 
 * Copyright (C) 2010  Pascal Fricke
 * Copyright (C) 2000-2010  Kirk Marat, The University of Manitoba
 * 
 * 
 * This file is part of nmr-fid library. 
 * nmr-fid library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * nmr-fid library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with cuteNMR.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * This code is based upon the libraries released by Dr Kirk Marat from the University
 * of Manitoba and cuteNMR:
 *      ftp://davinci.chem.umanitoba.ca/pub/marat/SpinWorks/source_library/
 *      https://sourceforge.net/projects/cutenmr/
 */
package uk.ac.ebi.nmr.fid.io;

import org.apache.commons.lang3.ArrayUtils;
import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;
import uk.ac.ebi.nmr.fid.Spectrum;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Reader for Bruker's fid file. The code is based in cuteNMR's code.
 *
 * @author Luis F. de Figueiredo
 *
 */
public class CuteNMRFidReader implements FidReader {

    private Acqu acquisition;
    private Proc processing;
    private InputStream fidInput;
    private double zeroFrequency = 0;

    public CuteNMRFidReader(InputStream fidInputS, Acqu acquisition) {
        this.acquisition = acquisition;
        this.fidInput = fidInputS;
    }

    public CuteNMRFidReader(InputStream fidInputS, Acqu acquisition, Proc processing) {
        this.acquisition = acquisition;
        this.fidInput = fidInputS;
        this.processing = processing;
    }

    public Spectrum read() throws Exception {



        DataInputStream in = null;

        in = new DataInputStream(new BufferedInputStream(fidInput));
        byte[] result = new byte[4];
        int[] fidInt = null;
        List<Integer> data = new ArrayList<Integer>();

        int totalBytesRead = 0;
        System.out.println(acquisition.is32Bit());
        if (acquisition.is32Bit()) {
            fidInt = new int[acquisition.getAquiredPoints()];

            // read each set of 4 bytes in the fid to get the intensity values
            while (in.available() >= 4) {
                // byteOrder == 1 => Big-Endian else Little-Endian
                // this has to do with the order of the bytes
                data.add((acquisition.getByteOrder() == 1)
                        ? (in.readByte()) << 24
                        | (in.readByte() << 16)
                        | (in.readByte() << 8)
                        | (in.readByte())
                        : (in.readByte())
                        | (in.readByte() << 8)
                        | (in.readByte() << 16)
                        | (in.readByte() << 24));

            }

            System.out.println(data.isEmpty());
            System.out.println("bytes left: " + in.available());
            in.close();


            //split imaginary from real part
//                List<Integer>realPoints = new ArrayList<Integer>();
//                List<Integer> imagPoints = new ArrayList<Integer>();
//                for( int i = 0 ; i < points.size()/2 ; i++ ){
//                    realPoints.add(points.get(i*2)); //even index numbers contain real fid
//                    imagPoints.add(points.get(i*2+1));//odd index numbers contain imaginary fid
//                }

        } else {
            //TODO read 64bit FIDs... and store it as double?
//            fidDouble = new double[acquisition.getAquiredPoints()];

        }

        System.out.println("Number of aquiredPoints and number of points in the FID:" + acquisition.getAquiredPoints() + " == "
                + data.size());

        System.out.println("Num bytes read: " + totalBytesRead);

        fidInt = ArrayUtils.toPrimitive(data.toArray(new Integer[data.size()]));
        double[] fid = new double[fidInt.length];
        for(int i =0; i<fidInt.length;i++)
            fid[i]=(double) fidInt[i];

        return new Spectrum(fid, acquisition);

    }
}
