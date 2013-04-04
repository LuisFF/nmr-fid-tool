package uk.ac.ebi.nmr.fid.io;

import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.nmr.fid.Fid;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 02/04/2013
 * Time: 16:02
 * To change this template use File | Settings | File Templates.
 */
public class ConnjurFidReaderTest {

    @Test
    public void testRead() throws Exception {
       ConnjurFidReader fidReader = new ConnjurFidReader(new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/" +
               "resources/examples/file_formats/bmse000109/1H/"));

        Fid fid = fidReader.read();
        Assert.assertEquals("Wrong number of data points", 19478, fid.getData().length);
        System.out.println("Nb of points: "+fid.getData().length);
    }
}
