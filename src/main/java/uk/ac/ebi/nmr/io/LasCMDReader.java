package uk.ac.ebi.nmr.io;

import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;

import javax.vecmath.Point3d;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 26/09/2012
 * Time: 17:58
 * To change this template use File | Settings | File Templates.
 */
public class LasCMDReader {
    BufferedReader bufferedReader;
    private final static Pattern REGEX_BONDS = Pattern.compile("(\\d+)\\s+(\\d+)"); //pair of atom index and bond order (i=7,8,...) - it may not work
    private final static Pattern REGEX_ATOM_TABLE = Pattern.compile("\\s?(\\w+)" + //element symbol (i=1)
            "\\s+[a-z0-9]+\\s+\\w+\\s+\\w+\\s+\\d+\\s+"+ // junk...
            "(-?\\d+\\.\\d+)\\s+\\d+\\s+\\d+\\s+"+ // charge (i=2)
            "(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)\\s+"+ // coordinates (i=3-5)
            "(\\d+)\\s+", Pattern.CASE_INSENSITIVE); // number of bonds
    private final static Pattern REGEX_COUNTS = Pattern.compile("(\\d+)\\s+(\\d+)\\s+(\\d+)"); //atom count
    private final static Pattern REGEX_COUNTS_TRIGGER = Pattern.compile("#Geometry, nAtoms, charge, spin multiplicity"); //Count trigger expression

    public LasCMDReader(Reader fileReader) {
        this.bufferedReader=new BufferedReader(fileReader);
    }

    public LasCMDReader(InputStream inputStream) {
        this(new InputStreamReader(inputStream));
    }

    public void setBufferedReader(InputStream inputReader) {
        setBufferedReader(new InputStreamReader(inputReader));
    }

    public void setBufferedReader(Reader fileReader) {
        this.bufferedReader = new BufferedReader(fileReader);
    }

    public boolean accepts(Class classObject) {
        Class[] interfaces = classObject.getInterfaces();
        for (int i = 0; i < interfaces.length; i++) {
            if (IChemFile.class.equals(interfaces[i])) return true;
            if (IChemSequence.class.equals(interfaces[i])) return true;
        }
        return false;
    }

    public void close() throws IOException {
        bufferedReader.close();
    }

    public <T extends IChemObject> T read(T object) throws IOException {
        if (object instanceof IChemSequence) {
            return (T) readChemSequence((IChemSequence) object);
        } else if (object instanceof IChemFile) {
            return (T) readChemFile((IChemFile) object);
        } else {
            throw new IOException("Object " + object.getClass().getName() + " is not supported");
        }
    }

    private IChemFile readChemFile(IChemFile chemFile) throws IOException {
        IChemSequence sequence = readChemSequence(chemFile.getBuilder().newInstance(IChemSequence.class));
        chemFile.addChemSequence(sequence);
        return chemFile;
    }

    private IChemSequence readChemSequence(IChemSequence sequence) throws IOException {
        IChemModel model = null;
        Matcher matcher;
        // Define the atomcontainer given as input in the Z-matrix
        model = sequence.getBuilder().newInstance(IChemModel.class);
        IAtomContainer container = model.getBuilder().newInstance(IAtomContainer.class);
        try {
            String line = bufferedReader.readLine();
            // store bond information
            int [][] adjacencyMatrix =null;

            // buffer ready does not read the last line
//            while (bufferedReader.ready() && (line != null)) {
            while (line != null) {
                if (REGEX_ATOM_TABLE.matcher(line).find()){
//                    System.out.println(line);
                    matcher=REGEX_ATOM_TABLE.matcher(line);
                    matcher.find();
//                    System.out.println(matcher.groupCount());
                    IAtom atom = container.getBuilder().newInstance(IAtom.class, matcher.group(1));
//                    atom.setCharge(Double.parseDouble(matcher.group(2)));
                    atom.setPoint3d(new Point3d(Double.parseDouble(matcher.group(3)),
                            Double.parseDouble(matcher.group(4)),
                            Double.parseDouble(matcher.group(5))));
                    container.addAtom(atom);
                    int nbOfBonds = Integer.parseInt(matcher.group(6));
                    int regionEnd = matcher.regionEnd();
                    int end = matcher.end();
                    // read bond information such as atom index and order
                    // store bond information because atoms might not exist yet...
                    // moreover, take care with duplicate information
                    for(int i = 0 ; i < nbOfBonds; i++)
                        if(REGEX_BONDS.matcher(line).region(end, regionEnd).find()){
                            matcher=REGEX_BONDS.matcher(line);
                            matcher.region(end, regionEnd).find();
                            end = matcher.end();
                            adjacencyMatrix[container.getAtomCount()-1][Integer.parseInt(matcher.group(1))-1]=
                                    Integer.parseInt(matcher.group(2));
//                            System.out.println(matcher.group(1) +" with order "+ matcher.group(2));

                        }
                }
                if(REGEX_COUNTS_TRIGGER.matcher(line).find()){
                    // TODO find a nicer way to extract the number of atoms
                    line =           bufferedReader.readLine();
                    line =           bufferedReader.readLine();

//                    System.out.println(line);
                    matcher=REGEX_COUNTS.matcher(line);
                    matcher.find();
                    adjacencyMatrix = new int[Integer.parseInt(matcher.group(1))][Integer.parseInt(matcher.group(1))];
                }

                line=bufferedReader.readLine();

            }
            IBond.Order[] orders = IBond.Order.values();
            // add bond information to the atom container
            // just needs the lower diagonal
            for (int i =1 ; i < adjacencyMatrix.length; i++)
                for (int j = 0; j < i; j++)
                    if (adjacencyMatrix[i][j]>0)
                        container.addBond(i, j, orders[adjacencyMatrix[i][j] - 1]);

            CDKHueckelAromaticityDetector.detectAromaticity(container);


        } catch (IOException exception) {
            throw new IOException("Error while reading general structure: " + exception.toString(), exception);
        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        // define the model as a molecule set
        // in order to store the different geometries
        IMoleculeSet moleculeSet = model.getBuilder().newInstance(IMoleculeSet.class);
        moleculeSet.addMolecule(model.getBuilder().newInstance(IMolecule.class, container));
        model.setMoleculeSet(moleculeSet);
        sequence.addChemModel(model);
        return sequence;
    }


}
