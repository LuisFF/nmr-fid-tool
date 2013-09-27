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

import edu.uchc.connjur.spectrumtranslator.UnexpectedDataException;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;
import uk.ac.ebi.nmr.fid.Spectrum;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;

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

    static double [] fid;
    static double [] pdata1r;
    static double [] pdata1i;

    @BeforeClass
    //TODO make this class available to all the tests of simply store the object
    public static void loadExternalData () {
        long startTime = System.currentTimeMillis();
        ObjectInputStream ois = null;
        try {
            ois = new ObjectInputStream(FidReader.class.getClassLoader()
                        .getResourceAsStream("data/bmse000109/1h-fid-raw.ser"));

            fid = (double []) ois.readObject();
            long endTime   = System.currentTimeMillis();
            System.out.println("Size of the fid:\t\t\t"+fid.length);
            System.out.println("Time reading fid object:\t\t"+(endTime - startTime)+" ms");
            startTime = System.currentTimeMillis();
            ois = new ObjectInputStream(FidReader.class.getClassLoader()
                    .getResourceAsStream("data/bmse000109/1h-pdata-1r.ser"));
            pdata1r = (double []) ois.readObject();
            endTime   = System.currentTimeMillis();
            System.out.println("Size of the 1r:\t\t\t"+pdata1r.length);
            System.out.println("Time reading 1r object:\t\t"+(endTime - startTime)+" ms");
            ois = new ObjectInputStream(FidReader.class.getClassLoader()
                    .getResourceAsStream("data/bmse000109/1h-pdata-1i.ser"));
            pdata1i = (double []) ois.readObject();
            endTime   = System.currentTimeMillis();
            System.out.println("Size of the 1i:\t\t\t"+pdata1i.length);
            System.out.println("Time reading 1i object:\t\t"+(endTime - startTime)+" ms");


        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }

    /**
     * Test Bruker DISP acquisition, 32 bit.
     * @throws IOException
     * @throws UnexpectedDataException
     */
    @Test
    public void testSimple1DFidReader() throws Exception {

//        File fidFile = new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/"+
//                "resources/examples/file_formats/bmse000109/1H/fid");

            Acqu acquisition = new BrukerAcquReader(this.getClass().
                    getClassLoader().getResourceAsStream("bmrb/1H/acqus")).read();

            FidReader fidReader = new Simple1DFidReader(new FileInputStream(this.getClass()
                    .getClassLoader().getResource("bmrb/1H/fid").getPath()),acquisition);
            Spectrum spepSpectrum = fidReader.read();
            Assert.assertNotNull("fid was not properly read", spepSpectrum);
            Assert.assertArrayEquals("fid was not properly read", fid, spepSpectrum.getFid(),1E-12);

    }
    /**
     * Test Bruker DISP acquisition, 32 bit. Without carrying about the FID
     * @throws IOException
     * @throws UnexpectedDataException
     */
    @Test
    public void testSimplePdataReader() throws Exception {

//        File fidFile = new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/"+
//                "resources/examples/file_formats/bmse000109/1H/fid");

        Acqu acquisition = new BrukerAcquReader(this.getClass().
                getClassLoader().getResourceAsStream("bmrb/1H/acqus")).read();
        Proc processing = new BrukerProcReader(this.getClass().
                getClassLoader().getResourceAsStream("bmrb/1H/pdata/1/procs"),acquisition).read();

        FidReader fidReader = new SimplePdataReader(new FileInputStream(this.getClass()
                .getClassLoader().getResource("bmrb/1H/pdata/1/1r").getPath()),new FileInputStream(this.getClass()
                .getClassLoader().getResource("bmrb/1H/pdata/1/1i").getPath()), acquisition,processing);
        Spectrum spepSpectrum = fidReader.read();
        Assert.assertNotNull("fid was not properly read", spepSpectrum);
        Assert.assertArrayEquals("fid was not properly read", pdata1r, spepSpectrum.getRealChannelData(),1E-12);
        Assert.assertArrayEquals("fid was not properly read", pdata1i, spepSpectrum.getImaginaryChannelData(),1E-12);

    }
    /**
     * Test Bruker DISP acquisition, 32 bit.
     * @throws IOException
     * @throws UnexpectedDataException
     */
    @Test
    public void testConnjurFidReader () throws Exception {
        Acqu acquisition = new BrukerAcquReader(this.getClass().
                getClassLoader().getResourceAsStream("bmrb/1H/acqu")).read();

        FidReader fidReader = new ConnjurFidReader(new File(this.getClass()
                .getClassLoader().getResource("bmrb/1H/").getFile()),acquisition);

        Spectrum spepSpectrum = fidReader.read();
        Assert.assertNotNull("fid was not properly read", spepSpectrum);
        Assert.assertArrayEquals("fid was not properly read", fid, spepSpectrum.getFid(),1E-12);

    }
}
