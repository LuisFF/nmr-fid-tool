package uk.ac.ebi.nmr.fid.io;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import uk.ac.ebi.nmr.fid.Acqu;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 02/04/2013
 * Time: 16:35
 * To change this template use File | Settings | File Templates.
 */
public class ConnjurAcquReaderTest {

    @Test
    public void testRead() throws Exception {
         AcquReader acquReader = new ConnjurAcquReader(new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/" +
                 "resources/examples/file_formats/bmse000109/1H/"));
        Acqu acqu = acquReader.read();

        Assert.assertNotNull("Acqu object is null",acqu);

        Assert.assertEquals("Wrong number of points", 19478, acqu.getAquiredPoints());
        System.out.println(acqu.getSpectralFrequency()); //        499.84
        System.out.println(acqu.getSpectralWidth());//             12.991109
        System.out.println(acqu.getDspGroupDelay());//             67.9855957


    }
}
