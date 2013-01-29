package uk.ac.ebi.nmr.fid.io;

import uk.ac.ebi.nmr.fid.Acqu;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 14/01/2013
 * Time: 14:12
 * To change this template use File | Settings | File Templates.
 */
public class AcquReader {

    private static File acquFile;


    // parameters from acqu
    private final static Pattern REGEXP_SFO1 = Pattern.compile("\\#\\#\\$SFO1= (-?\\d+\\.\\d+)"); //transmitter frequency
    private final static Pattern REGEXP_SFO2 = Pattern.compile("\\#\\#\\$SFO2= (-?\\d+\\.\\d+)"); //decoupler frequency
    private final static Pattern REGEXP_SFO3 = Pattern.compile("\\#\\#\\$SFO3= (\\d+\\.\\d+)"); //second decoupler frequency
    private final static Pattern REGEXP_O1 = Pattern.compile("\\#\\#\\$O1= (\\d+\\.\\d+)"); //frequency offset in Hz
    private final static Pattern REGEXP_SW = Pattern.compile("\\#\\#\\$SW= (\\d+\\.\\d+)"); //spectral width (ppm)
    private final static Pattern REGEXP_TD = Pattern.compile("\\#\\#\\$TD= (\\d+)"); //acquired points (real+imaginary)
    private final static Pattern REGEXP_DECIM = Pattern.compile("\\#\\#\\$DECIM= (-?\\d+)"); //DSP decimation factor
    private final static Pattern REGEXP_DSPFVS = Pattern.compile("\\#\\#\\$DSPFVS= (-?\\d+)"); //DSP firmware version
    private final static Pattern REGEXP_GRPDLY = Pattern.compile("\\#\\#\\$GRPDLY= (-?\\d+)"); //DSP group delay
    private final static Pattern REGEXP_BYTORDA = Pattern.compile("\\#\\#\\$BYTORDA= (\\d+)"); //byte order
    // variables not yet defined in Experiment
    private final static Pattern REGEXP_AQ_MODE = Pattern.compile("\\#\\#\\$AQ\\_mod= (\\d+)"); //acquisition mode
    private final static Pattern REGEXP_DIGMOD = Pattern.compile("\\#\\#\\$DIGMOD= (\\d+)"); //filter type
    private final static Pattern REGEXP_NS = Pattern.compile("\\#\\#\\$NS= (\\d+)"); //number of scans
    //TODO review REGEXP_PULPROG
    // examples of REGEXP_PULPROG : <zg> <cosydfph> <bs_hsqcetgpsi>; basically a word between < >
    private final static Pattern REGEXP_PULPROG = Pattern.compile("\\#\\#\\$PULPROG= (.+)"); //pulse program
    //TODO review REGEXP_NUC1
    // examples of REGEXP_NUC1 : <1H>; basically <isotope number + element>
    private final static Pattern REGEXP_NUC1 = Pattern.compile("\\#\\#\\$NUC1= (.+)"); // observed nucleus
    //TODO review REGEXP_INSTRUM
    // examples of REGEXP_INSTRUM : <amx500> ; basically <machine name>
    private final static Pattern REGEXP_INSTRUM = Pattern.compile("\\#\\#\\$INSTRUM= (.+)"); // instrument name
    private final static Pattern REGEXP_DTYPA = Pattern.compile("\\#\\#\\$DTYPA= (\\d+)"); //data type (0 -> 32 bit int, 1 -> 64 bit double)
    //TODO review REGEXP_SOLVENT
    // examples of REGEXP_SOLVENT : <DMSO> ; basically <solvent name>
    private final static Pattern REGEXP_SOLVENT = Pattern.compile("\\#\\#\\$SOLVENT= (.+)"); // instrument name
    //TODO review REGEXP_PROBHD
    // examples of REGEXP_PROBHD : <32> <>; basically <digit?>
    private final static Pattern REGEXP_PROBHD = Pattern.compile("\\#\\#\\$PROBHD= (.+)"); // probehead
    //TODO review REGEXP_ORIGIN
    // examples of REGEXP_ORIGIN : Bruker Analytik GmbH; basically a name
    private final static Pattern REGEXP_ORIGIN = Pattern.compile("\\#\\#\\$ORIGIN= (.+)"); // origin
    //TODO review REGEXP_OWNER
    // examples of REGEXP_OWNER : guest; basically the used ID
    private final static Pattern REGEXP_OWNER = Pattern.compile("\\#\\#\\$OWNER= (.+)"); // owner

    public AcquReader() {
    }

    public AcquReader(File acquFile) {
        this.acquFile = acquFile;
    }

    public AcquReader(String filename) throws FileNotFoundException {
        this.acquFile = new File(filename);
        // required parameters so far...
        // AquiredPoints: FidReader
        // SpectraWidth: FidReader
        // transmiterFreq: FidReader
        //
    }

    public Acqu read() throws Exception {
        Matcher matcher;
        Acqu acquisition = new Acqu();
        BufferedReader input = new BufferedReader(new FileReader(acquFile));
        String line = input.readLine();
        while (input.ready() && (line != null)) {
            if (REGEXP_SFO1.matcher(line).find()) {
                matcher = REGEXP_SFO1.matcher(line);
                matcher.find();
                acquisition.setTransmiterFreq(Double.parseDouble(matcher.group(1)));
            }
            if (REGEXP_SFO2.matcher(line).find()) {
                matcher = REGEXP_SFO2.matcher(line);
                matcher.find();
                acquisition.setDecoupler1Freq(Double.parseDouble(matcher.group(1)));
            }
            if (REGEXP_SFO3.matcher(line).find()) {
                matcher = REGEXP_SFO3.matcher(line);
                matcher.find();
                acquisition.setDecoupler2Feq(Double.parseDouble(matcher.group(1)));
            }
            if (REGEXP_O1.matcher(line).find()) {
                matcher = REGEXP_O1.matcher(line);
                matcher.find();
                acquisition.setFreqOffset(Double.parseDouble(matcher.group(1)));
            }
            if (REGEXP_SW.matcher(line).find()) {
                matcher = REGEXP_SW.matcher(line);
                matcher.find();
                acquisition.setSpectralWidth(Double.parseDouble(matcher.group(1)));
            }
            if (REGEXP_TD.matcher(line).find()) {
                matcher = REGEXP_TD.matcher(line);
                matcher.find();
                acquisition.setAquiredPoints(Integer.parseInt(matcher.group(1)));
            }
            if (REGEXP_DECIM.matcher(line).find()) {
                matcher = REGEXP_DECIM.matcher(line);
                matcher.find();
                acquisition.setDspDecimation(Integer.parseInt(matcher.group(1)));
            }
            if (REGEXP_DSPFVS.matcher(line).find()) {
                matcher = REGEXP_DSPFVS.matcher(line);
                matcher.find();
                acquisition.setDspFirmware(Integer.parseInt(matcher.group(1)));
            }
            if (REGEXP_GRPDLY.matcher(line).find()) {
                matcher = REGEXP_GRPDLY.matcher(line);
                matcher.find();
                acquisition.setDspGroupDelay(Integer.parseInt(matcher.group(1)));
            }
            if (REGEXP_BYTORDA.matcher(line).find()) {
                matcher = REGEXP_BYTORDA.matcher(line);
                matcher.find();
                acquisition.setByteOrder(Integer.parseInt(matcher.group(1)));
            }
            if (REGEXP_AQ_MODE.matcher(line).find()) {
                matcher = REGEXP_AQ_MODE.matcher(line);
                matcher.find();
                acquisition.setAcquisitionMode(Integer.parseInt(matcher.group(1)));
            }
            if (REGEXP_DIGMOD.matcher(line).find()) {
                matcher = REGEXP_DIGMOD.matcher(line);
                matcher.find();
                acquisition.setFilterType(Integer.parseInt(matcher.group(1)));
            }
            if (REGEXP_NS.matcher(line).find()) {
                matcher = REGEXP_NS.matcher(line);
                matcher.find();
                acquisition.setNumberOfScans(Integer.parseInt(matcher.group(1)));
            }
            if (REGEXP_PULPROG.matcher(line).find()) {
                matcher = REGEXP_PULPROG.matcher(line);
                matcher.find();
                acquisition.setPulseProgram(matcher.group(1));
            }
            if (REGEXP_NUC1.matcher(line).find()) {
                matcher = REGEXP_NUC1.matcher(line);
                matcher.find();
                acquisition.setObservedNucleus(matcher.group(1));
            }
            if (REGEXP_INSTRUM.matcher(line).find()) {
                matcher = REGEXP_INSTRUM.matcher(line);
                matcher.find();
                acquisition.setInstrumentName(matcher.group(1));
            }
            if (REGEXP_DTYPA.matcher(line).find()) {
                matcher = REGEXP_DTYPA.matcher(line);
                matcher.find();
                acquisition.set32Bit((Pattern.compile("0").matcher(matcher.group(1)).find()));
            }
            if (REGEXP_SOLVENT.matcher(line).find()) {
                matcher = REGEXP_SOLVENT.matcher(line);
                matcher.find();
                acquisition.setSolvent(matcher.group(1));
            }
            if (REGEXP_PROBHD.matcher(line).find()) {
                matcher = REGEXP_PROBHD.matcher(line);
                matcher.find();
                acquisition.setProbehead(matcher.group(1));
            }
            if (REGEXP_ORIGIN.matcher(line).find()) {
                matcher = REGEXP_ORIGIN.matcher(line);
                matcher.find();
                acquisition.setOrigin(matcher.group(1));
            }
            if (REGEXP_OWNER.matcher(line).find()) {
                matcher = REGEXP_OWNER.matcher(line);
                matcher.find();
                acquisition.setOwner(matcher.group(1));
            }
            line = input.readLine();
        }
        input.close();
        return acquisition;
    }


}
