package com.example.root.ffttest3;

import static com.example.root.ffttest3.Constants.LOG;

import android.util.Log;

import java.util.Arrays;

public class SymbolGeneration {
    public static short[] generate(short[] bits, int[] valid_carrier,
                                   int symreps, boolean preamble, Constants.SignalType sigType) {
        Log.e(LOG, "SymbolGeneration_generate");
        int numrounds = 0;
        if (valid_carrier.length > 0) {
            numrounds = (int) Math.ceil((double)bits.length/valid_carrier.length);
        }
        Log.e("sym",bits.length+","+valid_carrier.length+","+numrounds);

        int symlen = (Constants.Ns+Constants.Cp)*symreps + Constants.Gi;
//        int symlen = (Constants.Ns)*symreps + Constants.Gi;

        int siglen = symlen*numrounds;
        if (preamble) {
            siglen += ((Constants.preambleTime/1000.0)*Constants.fs)+Constants.ChirpGap;
        }
        short[] txsig = new short[siglen];

        int counter = 0;
        if (preamble) {
            // add preamble
            short[] preamble_sig = ChirpGen.preamble_s();
            for (Short s : preamble_sig) {
                txsig[counter++] = s;
            }
            counter += Constants.ChirpGap;
        }

        short[][] bit_list = new short[numrounds][valid_carrier.length];
        int bit_counter = 0;
        for (int i = 0; i < numrounds; i++) {
            int endpoint = bit_counter + valid_carrier.length-1;
            if (bit_counter + valid_carrier.length-1 > bits.length-1) {
                endpoint = bits.length-1;
            }
            short[] bits_seg = Utils.segment(bits,bit_counter,endpoint);
            if (Constants.DIFFERENTIAL && i > 0) {
                bits_seg = transform_bits(bit_list[i-1], bits_seg);
            }
            bit_list[i] = bits_seg;
            short[] symbol = generate_helper(
                    bits_seg,
                    valid_carrier,
                    symreps,
                    sigType);
            bit_counter += valid_carrier.length;

            for (Short s : symbol) {
                txsig[counter++] = s;
            }
        }
        Log.e("time",valid_carrier.length+","+bits.length+","+numrounds+","+(txsig.length/48000.0)+"");
        return txsig;
    }

    public static short[] generate2(short[] bits, int[] valid_carrier,
                                   int symreps, boolean preamble, Constants.SignalType sigType,
                                    int m_attempt) {
        Log.e(LOG, "SymbolGeneration_generate");
        int numrounds = 0;
        if (valid_carrier.length > 0) {
            numrounds = (int) Math.ceil((double)bits.length/valid_carrier.length);
        }
        Log.e("sym",bits.length+","+valid_carrier.length+","+numrounds+","+Constants.Ns+","+Constants.subcarrier_number_default);

        int symlen = (Constants.Ns+Constants.Cp)*symreps + Constants.Gi;

        int siglen = symlen*(numrounds+1);
        if (preamble) {
            if (Constants.DATA_NAISER) {
                siglen += ((Constants.preambleTime / 1000.0) * Constants.fs) + Constants.ChirpGap;
            }
            else {
                siglen += ((Constants.chirpPreambleTime/1000.0)*Constants.fs)+Constants.ChirpGap;
            }
        }
        short[] txsig = new short[siglen];

        int counter = 0;
        if (preamble) {
            // add preamble
            short[] preamble_sig = null;
            if (Constants.DATA_NAISER) {
                preamble_sig = ChirpGen.preamble_s();
            }
            else {
                preamble_sig = ChirpGen.chirp_get();
            }
            for (Short s : preamble_sig) {
                txsig[counter++] = s;
            }
            counter += Constants.ChirpGap;
        }

        // add training symbol
        short[][] bit_list = new short[numrounds+1][valid_carrier.length];

        short[] training_bits = null;
        if (Constants.subcarrier_number_default==20) {
            training_bits = Utils.segment(Constants.pn20_bits, 0, valid_carrier.length - 1);
        }
        else if (Constants.subcarrier_number_default==40) {
            training_bits = Utils.segment(Constants.pn40_bits, 0, valid_carrier.length - 1);
        }
        else if (Constants.subcarrier_number_default==60) {
            training_bits = Utils.segment(Constants.pn60_bits, 0, valid_carrier.length - 1);
        }
        else if (Constants.subcarrier_number_default==120) {
            training_bits = Utils.segment(Constants.pn120_bits, 0, valid_carrier.length - 1);
        }
        else if (Constants.subcarrier_number_default==300) {
            training_bits = Utils.segment(Constants.pn300_bits, 0, valid_carrier.length - 1);
        }
        else if (Constants.subcarrier_number_default==600) {
            training_bits = Utils.segment(Constants.pn600_bits, 0, valid_carrier.length - 1);
        }
        short[] symbol = generate_helper(
                training_bits,
                valid_carrier,
                symreps,
                sigType);
        for (Short s : symbol) {
            txsig[counter++] = s;
        }
        bit_list[0] = training_bits;

        int bit_counter = 0;
        Log.e("symbol", sigType.toString());
        Log.e("symbol", "# bits "+bits.length);
        Log.e("symbol", "# carriers "+valid_carrier.length);
        Log.e("symbol", "# symbols "+numrounds);

        String bitsWithPadding = "";
        String bitsWithoutPadding = "";
        String numberOfDataBits = "";

        for (int i = 0; i < numrounds; i++) {
            boolean oneMoreBin = i < bits.length%numrounds;

            int endpoint = (int)(bit_counter + Math.floor(bits.length/numrounds));
            if (!oneMoreBin) {
                endpoint -= 1;
            }

            short[] bits_seg = Utils.segment(bits,bit_counter,endpoint);
            numberOfDataBits += bits_seg.length+", ";
            bitsWithoutPadding += Utils.trim(Arrays.toString(bits_seg))+", ";

            short[] pad_bits = Utils.random_array(valid_carrier.length-bits_seg.length);
            Log.e("symbol", "sym "+i+": "+bits_seg.length+","+pad_bits.length);

            short[] tx_bits = Utils.concat_short(bits_seg,pad_bits);
            bitsWithPadding += Utils.trim(Arrays.toString(tx_bits))+", ";

            tx_bits = transform_bits(bit_list[i], tx_bits);
            for (int j = 0; j < tx_bits.length; j++) {
                bit_list[i+1][j] = tx_bits[j];
            }

            symbol = generate_helper(
                    tx_bits,
                    valid_carrier,
                    symreps,
                    sigType);
            bit_counter += bits_seg.length;

            for (Short s : symbol) {
                txsig[counter++] = s;
            }
        }
        int count = bitsWithPadding.length() - bitsWithPadding.replace(",", "").length();
        Log.e("time2",valid_carrier.length+","+bits.length+","+count+","+numrounds+","+(txsig.length/48000.0)+"");

        if (sigType.equals(Constants.SignalType.DataAdapt)) {
            FileOperations.writetofile(MainActivity.av, Utils.trim_end(bitsWithPadding),
                    Utils.genName(Constants.SignalType.BitsAdapt_Padding, m_attempt) + ".txt");
            FileOperations.writetofile(MainActivity.av, Utils.trim_end(bitsWithoutPadding),
                    Utils.genName(Constants.SignalType.BitsAdapt, m_attempt) + ".txt");
            FileOperations.writetofile(MainActivity.av, Utils.trim_end(numberOfDataBits),
                    Utils.genName(Constants.SignalType.Bit_Fill_Adapt, m_attempt) + ".txt");
        }
        else if (sigType.equals(Constants.SignalType.DataFull_1000_4000)) {
            FileOperations.writetofile(MainActivity.av, Utils.trim_end(bitsWithPadding),
                    Utils.genName(Constants.SignalType.BitsFull_1000_4000_Padding, m_attempt) + ".txt");
            FileOperations.writetofile(MainActivity.av, Utils.trim_end(bitsWithoutPadding),
                    Utils.genName(Constants.SignalType.BitsFull_1000_4000, m_attempt) + ".txt");
            FileOperations.writetofile(MainActivity.av, Utils.trim_end(numberOfDataBits),
                    Utils.genName(Constants.SignalType.Bit_Fill_1000_4000, m_attempt) + ".txt");
        }
        else if (sigType.equals(Constants.SignalType.DataFull_1000_2500)) {
            FileOperations.writetofile(MainActivity.av, Utils.trim_end(bitsWithPadding),
                    Utils.genName(Constants.SignalType.BitsFull_1000_2500_Padding, m_attempt) + ".txt");
            FileOperations.writetofile(MainActivity.av, Utils.trim_end(bitsWithoutPadding),
                    Utils.genName(Constants.SignalType.BitsFull_1000_2500, m_attempt) + ".txt");
            FileOperations.writetofile(MainActivity.av, Utils.trim_end(numberOfDataBits),
                    Utils.genName(Constants.SignalType.Bit_Fill_1000_2500, m_attempt) + ".txt");
        }
        else if (sigType.equals(Constants.SignalType.DataFull_1000_1500)) {
            FileOperations.writetofile(MainActivity.av, Utils.trim_end(bitsWithPadding),
                    Utils.genName(Constants.SignalType.BitsFull_1000_1500_Padding, m_attempt) + ".txt");
            FileOperations.writetofile(MainActivity.av, Utils.trim_end(bitsWithoutPadding),
                    Utils.genName(Constants.SignalType.BitsFull_1000_1500, m_attempt) + ".txt");
            FileOperations.writetofile(MainActivity.av, Utils.trim_end(numberOfDataBits),
                    Utils.genName(Constants.SignalType.Bit_Fill_1000_1500, m_attempt) + ".txt");
        }
        return txsig;
    }

    public static short[] transform_bits(short[] last_bits, short[] bits) {
        short[] newbits= new short[bits.length];
        for (int i = 0; i < bits.length; i++) {
            if (last_bits[i] != bits[i]) {
                newbits[i] = 1;
            }
        }
        return newbits;
    }

    public static short[] mod(int bound1, int bound2, int[] valid_carrier, short[] bits, int subnum, Constants.SignalType sigType) {
        double[][] mod = Modulation.pskmod(bits);
        if (bits.length < subnum) {
            bound2 = bound1+bits.length-1;
        }

        double[][] sig = new double[2][Constants.Ns];
        int counter=0;
        for (int i = bound1; i <= bound2; i++) {
            if (contains(valid_carrier, i)) {
                sig[0][i] = mod[0][counter];
                sig[1][i] = mod[1][counter++];
            }
        }

        double[][] symbol_complex = Utils.ifftnative2(sig);

        short[] symbol = new short[symbol_complex[0].length];
        double divval=(bound2-bound1)+1;
        if (sigType.equals(Constants.SignalType.Sounding)) {
            divval = divval/2;
        }
        for (int i = 0; i < symbol.length; i++) {
            symbol[i] = (short)((symbol_complex[0][i]/(double)divval)*32767.0);
        }

        return symbol;
    }

    // generate one symbol
    public static short[] generate_helper(short[] bits, int[] valid_carrier, int symreps, Constants.SignalType sigType) {
        int bound1=0;
        int bound2=0;
        int subnum=0;
        if (sigType.equals(Constants.SignalType.Sounding)) {
            bound1 = Constants.nbin1_default;
            bound2 = Constants.nbin2_default;
            subnum=Constants.subcarrier_number_chanest;
        }
        else if (sigType.equals(Constants.SignalType.DataAdapt)||
                sigType.equals(Constants.SignalType.DataFull_1000_4000)||
                sigType.equals(Constants.SignalType.DataFull_1000_2500)||
                sigType.equals(Constants.SignalType.DataFull_1000_1500)) {
            bound1 = valid_carrier[0];
            bound2 = valid_carrier[valid_carrier.length-1];
            subnum = bound2-bound1+1;
        }

        short[] symbol = mod(bound1, bound2, valid_carrier, bits, subnum, sigType);

        short[] flipped_symbol=null;
        flipped_symbol = mod(bound1, bound2, valid_carrier, Utils.flip(bits), subnum, sigType);

        int datacounter=0;
        short[] out = null;

        if (sigType.equals(Constants.SignalType.DataAdapt)||
            sigType.equals(Constants.SignalType.DataFull_1000_4000)||
            sigType.equals(Constants.SignalType.DataFull_1000_2500)||
            sigType.equals(Constants.SignalType.DataFull_1000_1500)) {
            int long_cp = Constants.Cp * symreps;
            short[] cp = new short[long_cp];
            for (int i = 0; i < long_cp; i++) {
                cp[i] = symbol[(symbol.length - long_cp - 1) + i];
            }

            out = new short[symbol.length*symreps + long_cp + Constants.Gi];

            for (int j = 0; j < cp.length; j++) {
                out[datacounter++] = cp[j];
            }
            for (int i = 0; i < symreps; i++) {
                for (int j = 0; j < symbol.length; j++) {
                    out[datacounter++] = symbol[j];
                }
            }
        }
        else if (sigType.equals(Constants.SignalType.Sounding)) {
            int long_cp = Constants.Cp;

            short[] cp_normal = new short[long_cp];
            for (int i = 0; i < long_cp; i++) {
                cp_normal[i] = symbol[(symbol.length - long_cp - 1) + i];
            }
            short[] cp_flipped = new short[long_cp];
            for (int i = 0; i < long_cp; i++) {
                cp_flipped[i] = flipped_symbol[(symbol.length - long_cp - 1) + i];
            }

            out = new short[(symbol.length+long_cp)*symreps + Constants.Gi];

            for (int i = 0; i < symreps; i++) {
                if (Utils.in(Constants.normal_syms,i)) {
                    for (int j = 0; j < cp_normal.length; j++) {
                        out[datacounter++] = cp_normal[j];
                    }
                    for (int j = 0; j < symbol.length; j++) {
                        out[datacounter++] = symbol[j];
                    }
                }
                else {
                    for (int j = 0; j < cp_normal.length; j++) {
                        out[datacounter++] = cp_flipped[j];
                    }
                    for (int j = 0; j < symbol.length; j++) {
                        out[datacounter++] = flipped_symbol[j];
                    }
                }
            }
        }

//        if (sigType.equals(Constants.SignalType.Sounding) && Constants.FLIP_SYMBOL) {
//            int scounter = 0;
//            for (int i = 0; i < symreps; i++) {
//                if (i == 1 || i == 2 || i == 3 || i == 4 || i==5) {
//                    for (int j = 0; j < symbol.length; j++) {
//                        out[scounter++] = symbol[j];
//                    }
//                }
//                else {
//                    for (int j = 0; j < symbol.length; j++) {
//                        out[scounter++] = flipped_symbol[j];
//                    }
//                }
//            }
//        }
//        else {
//            int counter = 0;
//        }

        return out;
    }

    public static boolean contains(int[] data, int k) {
        for (Integer i : data) {
            if (k==i) {
                return true;
            }
        }
        return false;
    }

    public static short[] rand_bits(int n) {
        Log.e(LOG, "SymbolGeneration_rand_bits");
        Constants.resetRandom();
        short[] arr = new short[n];
        for (int i = 0; i < arr.length; i++) {
            arr[i] = (short)Constants.random.nextInt(2);
        }
        return arr;
    }

    public static short[] getCodedBits(int n) {
        if (Constants.codeRate.equals(Constants.CodeRate.C1_2)) {
            return Utils.segment(Constants.data12,0,n-1);
        }
        else if (Constants.codeRate.equals(Constants.CodeRate.C2_3)) {
            return Utils.segment(Constants.data23,0,n-1);
        }
        return null;
    }
}
