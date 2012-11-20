package uk.ac.ebi.nmr.io;

import org.openscience.cdk.interfaces.IAtom;
import uk.ac.ebi.nmr.MagneticPropTool;
import uk.ac.ebi.nmr.SpinSystem;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Writer;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 04/10/2012
 * Time: 14:40
 * To change this template use File | Settings | File Templates.
 */
public class NMRDBSimulatorWriter {

    private static BufferedWriter bufferedWriter;
    public NMRDBSimulatorWriter(Writer writer) {

        this.bufferedWriter = new BufferedWriter(writer);

    }

    public NMRDBSimulatorWriter(OutputStream outputStream) {

        this(new OutputStreamWriter(outputStream));

    }

    public void write(SpinSystem spinSystem) throws IOException {
//        HashMap<IAtom, Integer> equivalentClasses = spinSystem.getChemEquivalenceMap();
        PrintWriter printWriter = new PrintWriter(bufferedWriter);

        for (IAtom atom : spinSystem.getAtomContainer().atoms()) {
            if(atom.getSymbol().equals("C")){
                for(IAtom connectedAtom : spinSystem.getAtomContainer().getConnectedAtomsList(atom)){
                    if(connectedAtom.getSymbol().equals("H")){
                        printWriter.printf("%d\t%.3f", spinSystem.getAtomContainer().getAtomNumber(connectedAtom),
                                (32.18375 - spinSystem.getMagneticTensors()
                                        [((Integer) atom.getProperty(MagneticPropTool.EQUIVALENT_CLASS)) - 1]));
                        for(IAtom atom2 : spinSystem.getAtomContainer().atoms()){
                            if(atom2.getSymbol().equals("C"))
                                if((Integer) atom.getProperty(MagneticPropTool.EQUIVALENT_CLASS) !=
                                        (Integer) atom2.getProperty(MagneticPropTool.EQUIVALENT_CLASS))
                                    for(IAtom connectedAtom2 : spinSystem.getAtomContainer().getConnectedAtomsList(atom2)){
                                        if(connectedAtom2.getSymbol().equals("H"))
                                            printWriter.printf("\t%d\t%.3f", spinSystem.getAtomContainer().getAtomNumber(connectedAtom2),
                                                    spinSystem.getConstantJMatrix()[
                                                            (Integer) atom.getProperty(MagneticPropTool.EQUIVALENT_CLASS)-1][
                                                            (Integer) atom2.getProperty(MagneticPropTool.EQUIVALENT_CLASS)-1]);
                                    }
                        }
                        printWriter.println();
                    }
                }
            }
        }
        printWriter.close();
    }

}
