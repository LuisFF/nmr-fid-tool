package uk.ac.ebi.nmr.fid.io;

import org.junit.Test;
import uk.ac.ebi.nmr.fid.io.BrukerReader;

import java.io.FileNotFoundException;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 21/09/2012
 * Time: 14:48
 * To change this template use File | Settings | File Templates.
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
