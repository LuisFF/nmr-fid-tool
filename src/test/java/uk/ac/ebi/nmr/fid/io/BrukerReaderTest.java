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

import org.junit.Test;

import java.io.FileNotFoundException;

/**
 * Test class for bruker reader class
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 21/09/2012
 * Time: 14:48
 *
 */
public class BrukerReaderTest {

    @Test
    public void testReadFID (){

        try {

            BrukerReader brukerReader  = new BrukerReader("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/"+
                    "resources/examples/file_formats/bmse000109/1H/fid");

//            brukerReader.read();


        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


    }

}
