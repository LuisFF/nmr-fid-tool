package uk.ac.ebi.nmr.io;

import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.rebond.RebondTool;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.io.DefaultChemObjectReader;
import org.openscience.cdk.io.formats.Gaussian03Format;
import org.openscience.cdk.io.formats.IResourceFormat;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.SaturationChecker;

import javax.vecmath.Point3d;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 25/09/2012
 * Time: 17:01
 * To change this template use File | Settings | File Templates.
 */
public class G09InputReader extends DefaultChemObjectReader {


    // Tables required to define the IAtomContainer in the input file
    private final static Pattern REGEX_ROOT = Pattern.compile("^#");
    private final static Pattern REGEX_EMPTY_LINE = Pattern.compile("^\\s*$");
    private final static Pattern REGEX_LINK= Pattern.compile("--link\\d+--",Pattern.CASE_INSENSITIVE);

    // assume that at least there is one digit in the left of the decimal place
    private final static Pattern REGEX_Z_MATRIX = Pattern
            .compile("\\s?(\\w+?)\\s+(-?\\d+\\.\\d*)\\s+(-?\\d+\\.\\d*)\\s+(-?\\d+\\.\\d*)");

    //TODO add regex for charge

    public static final String ROOT = "org.openscience.cdk.io.G09InputReader:root";
    public static final String HEADER = "org.openscience.cdk.io.G09InputReader:header";

    private BufferedReader input;
    private static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(G09InputReader.class);

    public G09InputReader(Reader reader) {
        input = new BufferedReader(reader);
    }

    public G09InputReader(InputStream input) {
        this(new InputStreamReader(input));
    }

    public G09InputReader() {
        this(new StringReader(""));
    }

    @TestMethod("testGetFormat")
    public IResourceFormat getFormat() {
        return Gaussian03Format.getInstance();
    }

    @TestMethod("testSetReader_Reader")
    public void setReader(Reader reader) throws CDKException {
        this.input = new BufferedReader(input);
    }

    @TestMethod("testSetReader_InputStream")
    public void setReader(InputStream input) throws CDKException {
        setReader(new InputStreamReader(input));
    }

    @TestMethod("testAccepts")
    public boolean accepts(Class classObject) {
        Class[] interfaces = classObject.getInterfaces();
        for (int i = 0; i < interfaces.length; i++) {
            if (IChemFile.class.equals(interfaces[i])) return true;
            if (IChemSequence.class.equals(interfaces[i])) return true;
        }
        return false;
    }

    public <T extends IChemObject> T read(T object) throws CDKException {
        if (object instanceof IChemSequence) {
            return (T) readChemSequence((IChemSequence) object);
        } else if (object instanceof IChemFile) {
            return (T) readChemFile((IChemFile) object);
        } else {
            throw new CDKException("Object " + object.getClass().getName() + " is not supported");
        }
    }

    @TestMethod("testClose")
    public void close() throws IOException {
        input.close();
    }

    private IChemFile readChemFile(IChemFile chemFile) throws CDKException {
        IChemSequence sequence = readChemSequence(chemFile.getBuilder().newInstance(IChemSequence.class));
        chemFile.addChemSequence(sequence);
        return chemFile;
    }

    /**
     * parse the file and read a chemicalSequence. At the moment, the information of the links is ignored.
     * @param sequence
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     */
    private IChemSequence readChemSequence(IChemSequence sequence) throws CDKException {
        IChemModel model = null;
        Matcher matcher;
        String header = null;
        String root = null;


        try {
            String line = input.readLine();
            // Define the atomcontainer given as input in the Z-matrix
            model = sequence.getBuilder().newInstance(IChemModel.class);
            IAtomContainer container = model.getBuilder().newInstance(IAtomContainer.class);

            // default values suggested by Eagon
            RebondTool rebonder = new RebondTool(2.0, 0.5, 0.5);
            boolean headerFlag = false;
            int emptyLines =0;

            while (input.ready() && (line != null)) {
                if(REGEX_EMPTY_LINE.matcher(line).find()){
                    emptyLines++;
                    line = input.readLine();
                }
                if (REGEX_Z_MATRIX.matcher(line).find()) {
                    //in previous versions of the reader, the symbol was defined based in the atomicNumber
                    matcher = REGEX_Z_MATRIX.matcher(line);
                    matcher.find();
                    IAtom atom = container.getBuilder().newInstance(IAtom.class, matcher.group(1));

                    atom.setPoint3d(new Point3d(Double.parseDouble(matcher.group(2)),
                            Double.parseDouble(matcher.group(3)),
                            Double.parseDouble(matcher.group(4))));
                    container.addAtom(atom);
                }
                // store header information as property in the atom container
                headerFlag |= REGEX_ROOT.matcher(line).find();
                if(headerFlag && emptyLines==1)
                    header = (header==null)?line: header+"\n"+line;

                // store root information as property in the atom container
                if(REGEX_ROOT.matcher(line).find())
                    root = (root==null)?line: root+"\n"+line;
                // do not store information for the links for the moment
                if (REGEX_LINK.matcher(line).find()){
                 break;
                }

                line = input.readLine();
            }
            // perform atom typing with AtomTypeFactory instead of AtomContainerManipulator
            // the latter doesn't work
            AtomTypeFactory factory = AtomTypeFactory.getInstance("org/openscience/cdk/config/data/jmol_atomtypes.txt"
                    , container.getBuilder());
            for (IAtom atom : container.atoms())
                factory.configure(atom);

            rebonder.rebond(container);

            CDKHueckelAromaticityDetector.detectAromaticity(container);
            new SaturationChecker().saturate(container);

            // define the model as a molecule set
            // in order to store the different geometries
            IMoleculeSet moleculeSet = model.getBuilder().newInstance(IMoleculeSet.class);
            moleculeSet.addMolecule(model.getBuilder().newInstance(IMolecule.class, container));
            model.setMoleculeSet(moleculeSet);

            model.setProperty(ROOT,root);
            model.setProperty(HEADER,header);

            sequence.addChemModel(model);
        } catch (IOException exception) {
            throw new CDKException("Error while reading general structure: " + exception.toString(), exception);
        }
        return sequence;
    }

}
