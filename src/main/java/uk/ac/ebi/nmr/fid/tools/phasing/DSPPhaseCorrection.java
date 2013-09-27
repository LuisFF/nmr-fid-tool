package uk.ac.ebi.nmr.fid.tools.phasing;

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Spectrum;

/**
 * Phase correction due to the DSP filter
 *
 * @author Luis F. de Figueiredo
 *
 */
public class DSPPhaseCorrection {


//    public DSPPhaseCorrection() {
//
//    }

    /**
     * Method to calculate the initial phase correction due to the digital filter. This method was based on the
     * information available in the dsp_phase of SpinWorks from Prof. Kirk Marat
     */
    public Spectrum dspPhaseCorrection(Spectrum spectrum) {
        double phase = 0;
        if (spectrum.getAcqu().getDspFirmware() == 10) {
            switch (spectrum.getAcqu().getDspDecimation()) {
                case 2:
                    phase = 16110.0;
                    break;
                case 3:
                    phase = 12060.0;
                    break;
                case 4:
                    phase = 23985.0;
                    break;
                case 6:
                    phase = 21270.0;
                    break;
                case 8:
                    phase = 24682.5;
                    break;
                case 12:
                    phase = 21735.0;
                    break;
                case 16:
                    phase = 25031.25;
                    break;
                case 24:
                    phase = 21967.5;
                    break;
                case 32:
                    phase = 25205.62;
                    break;
                case 48:
                    phase = 22083.75;
                    break;
                case 64:
                    phase = 25292.8;
                    break;
                case 96:
                    phase = 22141.87;
                    break;
                case 128:
                    phase = 25336.4;
                    break;
                case 192:
                    phase = 22170.93;
                    break;
                case 256:
                    phase = 25358.2;
                    break;
                case 384:
                    phase = 22185.46;
                    break;
                case 512:
                    phase = 25369.1;
                    break;
                case 768:
                    phase = 22192.73;
                    break;
                case 1024:
                    phase = 25374.55;
                    break;
                case 1536:
                    phase = 22196.36;
                    break;
                case 2048:
                    phase = 25377.27;
                    break;
                default:
                    phase = 0;
                    break;

            }
        } else if (spectrum.getAcqu().getDspFirmware() == 11) {
            switch (spectrum.getAcqu().getDspDecimation()) {
                case 2:
                    phase = 16560.0;
                    break;
                case 3:
                    phase = 13140.0;
                    break;
                case 4:
                    phase = 17280.0;
                    break;
                case 6:
                    phase = 18060.0;
                    break;
                case 8:
                    phase = 19170.0;
                    break;
                case 12:
                    phase = 25020.0;
                    break;
                case 16:
                    phase = 26010.0;
                    break;
                case 24:
                    phase = 25260.0;
                    break;
                case 32:
                    phase = 26190.0;
                    break;
                case 48:
                    phase = 25380.0;
                    break;
                case 64:
                    phase = 26280.0;
                    break;
                case 96:
                    phase = 25440.0;
                    break;
                case 128:
                    phase = 26100.0;
                    break;
                case 192:
                    phase = 25680.0;
                    break;
                case 256:
                    phase = 26010.0;
                    break;
                case 384:
                    phase = 25800.0;
                    break;
                case 512:
                    phase = 25965.0;
                    break;
                case 768:
                    phase = 25860.0;
                    break;
                case 1024:
                    phase = 25942.5;
                    break;
                case 1536:
                    phase = 25890.0;
                    break;
                case 2048:
                    phase = 25931.25;
                    break;
                default:
                    phase = 0;
                    break;

            }

        } else if (spectrum.getAcqu().getDspFirmware() == 12) {
            switch (spectrum.getAcqu().getDspDecimation()) {
                case 2:
                    phase = 16560.0;
                    break;
                case 3:
                    phase = 13140.0;
                    break;
                case 4:
                    phase = 17280.0;
                    break;
                case 6:
                    phase = 18060.0;
                    break;
                case 8:
                    phase = 19170.0;
                    break;
                case 12:
                    phase = 25020.0;
                    break;
                case 16:
                    phase = 25785.0;
                    break;
                case 24:
                    phase = 25260.0;
                    break;
                case 32:
                    phase = 25965.0;
                    break;
                case 48:
                    phase = 25380.0;
                    break;
                case 64:
                    phase = 26055.0;
                    break;
                case 96:
                    phase = 25440.0;
                    break;
                case 128:
                    phase = 26100.0;
                    break;
                case 192:
                    phase = 25680.0;
                    break;
                case 256:
                    phase = 26010.0;
                    break;
                case 384:
                    phase = 25800.0;
                    break;
                case 512:
                    phase = 25965.0;
                    break;
                case 768:
                    phase = 25860.0;
                    break;
                case 1024:
                    phase = 25942.5;
                    break;
                case 1536:
                    phase = 25890.0;
                    break;
                case 2048:
                    phase = 25931.25;
                    break;
                default:
                    phase = 0;
                    break;

            }
        } else if ((spectrum.getAcqu().getDspFirmware() >= 20 && spectrum.getAcqu().getDspFirmware() <= 23) ||
                spectrum.getAcqu().getAcquisitionMode().equals(Acqu.AcquisitionMode.CUSTOM_DISP)) {
            phase = spectrum.getAcqu().getDspGroupDelay() * 360;
        } else {
            phase = 0;
            System.err.println("[WARNING] It is not possible to identify the DSP used;");
        }
        double [] fid = new double[spectrum.getFid().length];

        for (int i = 0; i < spectrum.getFid().length/2; i++) {
            double phaseAngle = 2 * Math.PI / 360 * i / (spectrum.getFid().length/2) * phase;
            // real channel are even positions
            fid[i*2] = spectrum.getFid()[i*2] * Math.cos(phaseAngle) - spectrum.getFid()[i*2+1] * Math.sin(phaseAngle);
            // imaginary channel are off positions
            fid[i*2+1] = spectrum.getFid()[i*2] * Math.sin(phaseAngle) + spectrum.getFid()[i*2+1] * Math.cos(phaseAngle);
            
        }

        return new Spectrum(fid,spectrum.getAcqu(),spectrum.getProc());
    }
}