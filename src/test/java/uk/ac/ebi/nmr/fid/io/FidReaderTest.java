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

import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Spectrum;

import java.io.File;
import java.io.IOException;

/**
 * Test class for fid reader
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 14/01/2013
 * Time: 18:41
 *
 */
public class FidReaderTest {

    @Test
    public void testCuteNMRFidReader() throws IOException{
        File fidFile = new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/"+
                "resources/examples/file_formats/bmse000109/1H/fid");

        try {
            Acqu acquisition = new BrukerAcquReader("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/"+
                    "resources/examples/file_formats/bmse000109/1H/acqu").read();
            FidReader fidReader = new CuteNMRFidReader(this.getClass().getResourceAsStream(fidFile.getPath()),acquisition);
            Spectrum spepSpectrum = fidReader.read();
            Assert.assertNotNull("fid was not properly read" ,spepSpectrum);

        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


    }
}
