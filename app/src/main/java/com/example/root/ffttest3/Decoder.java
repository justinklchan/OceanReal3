package com.example.root.ffttest3;

import static com.example.root.ffttest3.Constants.LOG;
import static com.example.root.ffttest3.Constants.valid_carrier_preamble;

import android.content.Context;
import android.util.Log;

public class Decoder {

    public static void extract(Context cxt, double[] rec, int[] valid_bins) {
        Log.e(Constants.LOG, "SignalProcessing_extract");
        double[] preamble = ChirpGen.preamble_d();
        double[] filt = Utils.filter(rec);
//        int start_point = Utils.xcorr(preamble, filt, rec, rec.length);
        int start_point=0;
        start_point += preamble.length+Constants.ChirpGap+1;

        int rx_start = start_point;
        int rx_end = start_point+(Constants.sym_len*Constants.Nsyms)-1;
        int rx_len = (rx_end-rx_start)+1;
        Log.e(Constants.LOG, ">>"+rec.length+","+rx_start+","+rx_end+","+rx_len);
        if (rx_end-1 > rec.length || rx_start < 0) {
            Utils.log("Error extracting preamble from data signal");
            return;
        }
        double[] rxsig = Utils.segment(rec, rx_start, rx_end);

        short[] bits = SymbolGeneration.rand_bits(valid_bins.length*Constants.Nsyms);
        short[] txsig = SymbolGeneration.generate(bits, valid_bins,Constants.data_symreps, false,
                Constants.SignalType.DataAdapt);
        decode(rxsig, Utils.convert(txsig));
    }

    public static void test_decode(Context cxt) {
        long t1 = System.currentTimeMillis();

        double[] rx_file=FileOperations.readrawasset(cxt, R.raw.test_rx,30000);
        Log.e(LOG,"rx read " +(System.currentTimeMillis()-t1));

        t1 = System.currentTimeMillis();
        // tx file should not have warmup, preamble or first gap
        double[] tx_file=FileOperations.readrawasset(cxt, R.raw.test_tx,1);
        Log.e(LOG,"tx read " +(System.currentTimeMillis()-t1));

        decode(rx_file,tx_file);
    }

    public static void decode(double[] rx_file, double[] tx_file) {
        Log.e(Constants.LOG, "SignalProcessing_decode "+rx_file.length+","+tx_file.length);
        long t1 = System.currentTimeMillis();

        int rx_idx = 0;
        int tx_idx =  0;

        short[][] bits_rx = new short[Constants.Nsyms][valid_carrier_preamble.length];
        short[][] bits_tx = new short[Constants.Nsyms][valid_carrier_preamble.length];

        double sym_level_ber = 0;
//        for (int sym = 0; sym < Constants.Nsyms; sym++) {
        for (int sym = 0; sym < 24; sym++) {
            int len_rx = 0;
            if (Constants.eqMethod == Constants.EqMethod.Freq) {
                len_rx = Constants.Ns;
            }
            else if (Constants.eqMethod == Constants.EqMethod.Time) {
                len_rx = Constants.Ns + Constants.tap_num - 1;
            }

            int tx_start = tx_idx+Constants.Cp;
            int tx_end = tx_idx + Constants.Cp + Constants.Ns - 1;
            int rx_start = (rx_idx-Constants.sync_offset)+Constants.Cp;
            int rx_end = (rx_idx-Constants.sync_offset)+Constants.Cp+Constants.Ns-1;

            if (tx_end-1 > tx_file.length) {
                Utils.log("Error extracting tx symbol from decoder");
                break;
            }
            if (rx_end-1 > rx_file.length) {
                Utils.log("Error extracting rx symbol from decoder");
                break;
            }

            double[] symbol_tx = Utils.segment(tx_file, tx_start, tx_end);
            double[] symbol_rx = Utils.segment(rx_file, rx_start, rx_end);

//            if (Constants.pilots.contains(sym)) {
            if (sym==0) {
//                t1 = System.currentTimeMillis();
                Equalizer.equalizer_estimation(symbol_tx, symbol_rx);
//                Log.e(Constants.LOG,"train "+(System.currentTimeMillis()-t1));
            }

//            t1 = System.currentTimeMillis();
            double[][] spec_rx = null;
            double[] symbol_pred = null;
            if (Constants.eqMethod == Constants.EqMethod.Freq) {
                spec_rx = Equalizer.equalizer_recover_freq(symbol_rx, sym);
            }
            else if (Constants.eqMethod == Constants.EqMethod.Time) {
                //symbol_pred = Equalizer.equalizer_recover_time(symbol_rx);
                spec_rx = Utils.fftcomplexoutnative_double(symbol_pred, Constants.Ns);
            }

            double[][] spec_tx = Utils.fftcomplexoutnative_double(symbol_tx, Constants.Ns);

//            t1 = System.currentTimeMillis();
            bits_rx[sym] = Modulation.pskdemod(spec_rx, valid_carrier_preamble);
            bits_tx[sym] = Modulation.pskdemod(spec_tx, valid_carrier_preamble);

//            double ber = ber(bits_rx[sym],bits_tx[sym]);
//            sym_level_ber += ber;
//            Log.e(Constants.LOG,String.format("demod %d %.2f",sym,ber));

            tx_idx += Constants.sym_len;
            rx_idx += Constants.sym_len;
        }
        Log.e("timer2",""+(System.currentTimeMillis()-t1));

        Log.e(LOG, "BER by sym "+sym_level_ber);
        ber_by_freq(bits_rx, bits_tx);
    }

    public static void ber_by_freq(short[][] all_bits_rx, short[][] all_bits_tx) {
        Utils.log("SignalProcessing_ber_by_freq");
        double mean_err=0;
        for (int i = 0; i < valid_carrier_preamble.length; i++) {
            short[] bits_rx = new short[Constants.Nsyms];
            short[] bits_tx = new short[Constants.Nsyms];
            for (int sym = 0; sym < Constants.Nsyms; sym++) {
                bits_rx[sym] = all_bits_rx[sym][i];
                bits_tx[sym] = all_bits_tx[sym][i];
            }
            double err = ber(bits_rx,bits_tx);
            mean_err += err;
            Utils.log(String.format("%d %d %.2f", valid_carrier_preamble[i],Constants.f_seq.get(valid_carrier_preamble[i]),err));
        }
        Utils.log( "Mean BER "+(mean_err/ valid_carrier_preamble.length));
    }

    public static double ber(short[] bits_rx, short[] bits_tx) {
        int matches=0;
        for (int i = 0; i < bits_rx.length; i++) {
            matches += (bits_rx[i] == bits_tx[i]) ? 1 : 0;
        }
        double acc = (double)matches /  bits_rx.length;
        return 1-acc;
    }
}
