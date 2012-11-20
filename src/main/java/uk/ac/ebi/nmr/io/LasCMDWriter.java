package uk.ac.ebi.nmr.io;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 01/08/2012
 * Time: 09:27
 * To change this template use File | Settings | File Templates.
 */
public class LasCMDWriter {
    private IAtomContainer atomContainer = null;
    private FileWriter fileWriter = null;
    private BufferedWriter buffWriter = null;

    /*
    MMtypes used by Allouch
    Here are only the ones that are not equal to the element symbol or that have various variants
     */
    private static final String C="sp2 C carbonyl group \n";
    private static final String CA="sp2 C pure aromatic (benzene)\n";
    private static final String CB="sp2 aromatic C, 5&6 membered ring junction\n";
    private static final String CC="sp2 aromatic C, 5 memb. ring HIS\n";
    private static final String CD="sp2 C atom in the middle of: C=CD-CD=C\n";
    private static final String CK="sp2 C 5 memb.ring in purines\n";
    private static final String CM="sp2 C  pyrimidines in pos. 5 & 6\n";
    private static final String CN="sp2 C aromatic 5&6 memb.ring junct.(TRP)\n";
    private static final String CQ="sp2 C in 5 mem.ring of purines between 2 N\n";
    private static final String CT="sp3 aliphatic C\n";
    private static final String CV="sp2 arom. 5 memb.ring w/1 N and 1 H (HIS)\n";
    private static final String CW="sp2 arom. 5 memb.ring w/1 N-H and 1 H (HIS)\n";
    private static final String CSTAR="sp2 arom. 5 memb.ring w/1 subst. (TRP)\n";
    private static final String CY="nitrile C (Howard et al.JCC,16,243,1995)\n";
    private static final String CZ="sp C (Howard et al.JCC,16,243,1995)\n";
    private static final String C0="calcium\n";
    private static final String H="H bonded to nitrogen atoms\n";
    private static final String HC="H aliph. bond. to C without electrwd.group\n";
    private static final String H1="H aliph. bond. to C with 1 electrwd. group\n";
    private static final String H2="H aliph. bond. to C with 2 electrwd. group\n";
    private static final String H3="H aliph. bond. to C with 3 electrwd. group\n";
    private static final String HA="H arom. bond. to C without elctrwd. groups\n";
    private static final String H4="H arom. bond. to C with 1 electrwd. group\n";
    private static final String H5="H arom.at C with 2 elctrwd. gr,+HCOO group\n";
    private static final String HO="hydroxyl group\n";
    private static final String HS="hydrogen bonded to sulphur (pol?)\n";
    private static final String HW="H in TIP3P water\n";
    private static final String HP="H bonded to C next to positively charged gr\n";
    private static final String HZ="H bond sp C (Howard et al.JCC,16,243,1995)\n";
    private static final String IM="assumed to be Cl- (ion minus)\n";
    private static final String IB="big ion w/ waters' for vacuum (Na+, 6H2O)\n";
    private static final String MG="magnesium\n";
    private static final String N="sp2 nitrogen in amide groups\n";
    private static final String NA="sp2 N in 5 memb.ring w/H atom (HIS)\n";
    private static final String NB="sp2 N in 5 memb.ring w/LP (HIS,ADE,GUA)\n";
    private static final String NC="sp2 N in 6 memb.ring w/LP (ADE,GUA)\n";
    private static final String N2="sp2 N in amino groups\n";
    private static final String N3="sp3 N for charged amino groups (Lys, etc)\n";
    private static final String NT="sp3 N for amino groups amino groups \n";
    private static final String NSTAR="sp2 N \n";
    private static final String NY="nitrile N (Howard et al.JCC,16,243,1995)\n";
    private static final String O="carbonyl group oxygen\n";
    private static final String O2="carboxyl and phosphate group oxygen\n";
    private static final String OW="oxygen in TIP3P water\n";
    private static final String OH="oxygen in hydroxyl group\n";
    private static final String OS="ether and ester oxygen\n";
    private static final String S="S in disulfide linkage,pol:JPC,102,2399,98\n";
    private static final String SH="S in cystine\n";
    private static final String CU="copper\n";
    private static final String FE="iron\n";
    private static final String IP="assumed to be Na+ (ion plus)\n";
    private static final String Na="Na+, ions pol:J.PhysC,11,1541,(1978)\n";


    public LasCMDWriter(FileWriter fileWriter) {
        this.fileWriter = fileWriter;
        this.buffWriter = new BufferedWriter(fileWriter);
    }
    public LasCMDWriter(BufferedWriter buffWriter) {
        this.buffWriter = buffWriter;
    }


    public void write(IAtomContainer atomContainer) throws IOException {
        this.atomContainer = atomContainer;
        String out = "#RunType = Energy, Optimization, MD, MDConfo, REMDConfo\n" +
                "RunType=MDConfo\n" +
                "#Model = MM , Mopac , Orca or FireFly\n" +
                "Model=MM\n" +
                "SEKeys=PM6\n" +
                "#SEKeys=AM1\n" +
                "mopacCommand=/homes/ldpf/applications/mopac/MOPAC2012.exe\n" +
                "orcaCommand=orca\n" +
                "fireflyCommand=firefly\n" +
                "\n" +
                "#Confo\n" +
                "gaussianCommand=g03\n" +
                "fireflyCommand=firefly\n" +
                "numberOfGeometries=5000\n" +
                "tolEnergy=0.01\n" +
                "tolDistance=0.01\n" +
                "ConfoOptMM=TRUE\n" +
                "ConfoOptMopac=TRUE\n" +
                "ConfoOptMopacMethod=PM6 GNORM=0.001\n" +
                "ConfoOptFireFly=FALSE\n" +
                "# remove # if post processing required\n" +
                "#mopacKeywordsPost=PM6\n" +
                "gaussianKeywordsPost=B3LYP/6-31G*\n" +
                "#fireflyKeywordsPost=AM1\n" +
                "\n" +
                "#MM\n" +
                "# AMBER, UFF(not implemented), PAIRWISE\n" +
                "ForceFieldType=0\n" +
                "ForceFieldUseBond=TRUE\n" +
                "ForceFieldUseBend=TRUE\n" +
                "ForceFieldUseDihedral=TRUE\n" +
                "ForceFieldUseImproper=FALSE\n" +
                "ForceFieldUseNonBonded=TRUE\n" +
                "ForceFieldUseHydrogenBonded=FALSE\n" +
                "ForceFieldUsecoulomb=TRUE\n" +
                "ForceFieldUseVanderWals=TRUE\n" +
                "#  NOCONSTRAINTS = 0, BONDSCONSTRAINTS = 1, BONDSANGLESCONSTRAINTS = 2\n" +
                "ForceFieldConstraints=1\n" +
                "\n" +
                "#MD\n" +
                "updateFrequency=5\n" +
                "#Time in ps\n" +
                "heatTime = 0\n" +
                "equiTime = 0\n" +
                "runTime = 15\n" +
                "coolTime = 0\n" +
                "timeExchange = 0.01\n" +
                "heatTemp = 0\n" +
                "runTemp = 5000\n" +
                "runTempMax = 5000\n" +
                "nTemperatures = 10\n" +
                "#in fs\n" +
                "stepSize = 0.5\n" +
                "#  VERLET = 0, BEEMAN = 1, STOCHASTIC = 2\n" +
                "integrator = 0\n" +
                "#  NONE = 0, ANDERSEN = 1, BERENDSEN = 2, BUSSI = 3\n" +
                "thermostat = 0\n" +
                "friction=40\n" +
                "collide = 20\n" +
                "\n" +
                "#QuasiNewton\n" +
                "useQuasiNewton = TRUE\n" +
                "quasiNewtonMaxIterations = 20000\n" +
                "quasiNewtonUpdateFrequency = 100\n" +
                "quasiNewtonEpsilon  = 0.0001\n" +
                "quasiNewtonTolerence = 1e-16\n" +
                "quasiNewtonMaxLines =  25\n" +
                "\n" +
                "#ConjugateGradient\n" +
                "useConjugateGradient = FALSE\n" +
                "conjugateGradientGradientNorm = 1e-3\n" +
                "conjugateGradientMaxIterations = 100\n" +
                "conjugateGradientUpdateFrequency = 1\n" +
                "conjugateGradientMaxLines = 25\n" +
                "conjugateGradientInitialStep = 0.001\n" +
                "# 1 : Hestenes Stiefel,  2 : Fletcher Reeves, 3 : Polak Ribiere, 4 : Wolf Powell\n" +
                "conjugateGradientMethod = 1\n" +
                "\n" +
                "#Geometry, nAtoms, charge, spin multiplicity. For each atom : symbol, MMType, pdbType,"+
                " residueName, numResidue, charge, layer, x(Ang),y,z, nconn, num1, type1, num2, type2,...\n" +
                "Geometry\n"+
                atomContainer.getAtomCount()+" "+
               (int) AtomContainerManipulator.getTotalCharge(atomContainer)+" 1\n";

         //TODO check how to get the multiplicity from a molecule

        for(IAtom atom : atomContainer.atoms()){
            //symbol
            out += " "+atom.getSymbol()+" ";
            //MMType???
            out += atom.getSymbol()+" ";
            //pdbType???
            out += atom.getSymbol()+" ";
            //residueName???
            out += atom.getSymbol()+" ";
            //numResidue???
            out += "0 ";
            //charge
            out += (atom.getCharge()==null)? "0.00000 ": atom.getCharge()+" ";
            //layer????
            out += "2 0 ";
            // 3D coordinates
            out+=atom.getPoint3d().x+" "+atom.getPoint3d().y+" "+atom.getPoint3d().z+"  ";
            // number of neighbour atoms
            out+=atomContainer.getConnectedAtomsCount(atom)+" ";
            for (IAtom neighbourAtom : atomContainer.getConnectedAtomsList(atom)){
                out+=" "+(atomContainer.getAtomNumber(neighbourAtom)+1);
                out+= (atomContainer.getBond(neighbourAtom, atom).getOrder().equals(IBond.Order.SINGLE))? " 1":
                        (atomContainer.getBond(neighbourAtom, atom).getOrder().equals(IBond.Order.DOUBLE))? " 2":
                                (atomContainer.getBond(neighbourAtom, atom).getOrder().equals(IBond.Order.DOUBLE))? " 3":" 4";

            }
            out+="\n";
        }
        buffWriter.write(out);
        buffWriter.close();
    }
}

