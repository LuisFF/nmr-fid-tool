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

import java.io.File;

/**
 * Test class for connjur acqu reader
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 02/04/2013
 * Time: 16:35
 *
 */
public class ConnjurAcquReaderTest {

    @Test
    public void testRead() throws Exception {
         AcquReader acquReader = new ConnjurAcquReader(new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/" +
                 "resources/examples/file_formats/bmse000109/1H/"), Acqu.Spectrometer.BRUKER);
        Acqu acqu = acquReader.read();

        Assert.assertNotNull("Acqu object is null",acqu);

        Assert.assertEquals("Wrong number of points", 19478, acqu.getAquiredPoints());
        System.out.println(acqu.getSpectralFrequency()); //        499.84
        System.out.println(acqu.getSpectralWidth());//             12.991109
        System.out.println(acqu.getDspGroupDelay());//             67.9855957


    }
}
