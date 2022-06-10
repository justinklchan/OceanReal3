package com.example.root.ffttest3;

import android.Manifest;
import android.app.Activity;
import android.app.FragmentManager;
import android.content.ContentResolver;
import android.content.Context;
import android.content.Intent;
import android.content.SharedPreferences;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.media.AudioTrack;
import android.net.Uri;
import android.os.AsyncTask;
import android.os.Build;
import android.os.Bundle;
import android.preference.PreferenceManager;
import android.provider.Settings;
import android.text.Editable;
import android.text.TextWatcher;
import android.text.method.ScrollingMovementMethod;
import android.util.Log;
import android.view.View;
import android.view.Window;
import android.view.WindowManager;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.CompoundButton;
import android.widget.EditText;
import android.widget.Spinner;
import android.widget.Switch;
import android.widget.TextView;
import android.widget.Toast;

import net.yslibrary.android.keyboardvisibilityevent.KeyboardVisibilityEvent;
import net.yslibrary.android.keyboardvisibilityevent.KeyboardVisibilityEventListener;

import androidx.appcompat.app.AppCompatActivity;
import androidx.core.app.ActivityCompat;
import androidx.core.widget.NestedScrollView;

import com.jjoe64.graphview.GraphView;

import java.util.ArrayList;

public class MainActivity extends AppCompatActivity implements SensorEventListener {
    String[] perms = new String[]{Manifest.permission.RECORD_AUDIO, Manifest.permission.WRITE_EXTERNAL_STORAGE, Manifest.permission.READ_EXTERNAL_STORAGE,Manifest.permission.WRITE_SETTINGS};

    private static SensorManager sensorManager;
    private Sensor accelerometer;
    private Sensor gyroscope;
    static Activity av;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
        getWindow().setFlags(WindowManager.LayoutParams.FLAG_FULLSCREEN,
                WindowManager.LayoutParams.FLAG_FULLSCREEN);
        // Example of a call to a native method
        uiSetup();
        Constants.cResolver = getContentResolver();

        Constants.window = getWindow();

        sensorManager = (SensorManager) getSystemService(Context.SENSOR_SERVICE);
        accelerometer = sensorManager.getDefaultSensor(Sensor.TYPE_LINEAR_ACCELERATION);
        gyroscope = sensorManager.getDefaultSensor(Sensor.TYPE_GYROSCOPE);
        sensorManager.registerListener(this, accelerometer, SensorManager.SENSOR_DELAY_FASTEST);
        sensorManager.registerListener(this, gyroscope, SensorManager.SENSOR_DELAY_FASTEST);

        KeyboardVisibilityEvent.setEventListener(
                this,
                new KeyboardVisibilityEventListener() {
                    @Override
                    public void onVisibilityChanged(boolean isOpen) {
                        FullScreencall();
                    }
                });

        ActivityCompat.requestPermissions(this,
                perms,
                1234);
        Constants.setup(this);
//        double[] out=ChirpGen.preamble_d();
//        testme();
//        String out = Utils.convtest("0101");
//        Log.e("conv",out);
//        Tests.corr_test();
//        Tests.test_decode();
//        Tests.naiser_test();
//        Tests.sounding_test();
//        ChirpGen.sounding_signal_s();
//        Tests.bin_filling();
//        Tests.feedback_test();
//        checkSystemWritePermission();
//        Constants.setBrightness(255);
//        Tests.sounding_test2(this);
//        Tests.feedback_test();
    }

    private boolean checkSystemWritePermission() {
        boolean retVal = true;
        if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.M) {
            retVal = Settings.System.canWrite(this);
//            Log.d(TAG, "Can Write Settings: " + retVal);
            if(retVal){
//                Toast.makeText(this, "Write allowed :-)", Toast.LENGTH_LONG).show();
            }else{
//                Toast.makeText(this, "Write not allowed :-(", Toast.LENGTH_LONG).show();
//                FragmentManager fm = getFragmentManager();
//                PopupWritePermission dialogFragment = new PopupWritePermission();
//                dialogFragment.show(fm, getString(R.string.popup_writesettings_title));

                openAndroidPermissionsMenu();
            }
        }
        return retVal;
    }

    private void openAndroidPermissionsMenu() {
        Intent intent = new Intent(Settings.ACTION_MANAGE_WRITE_SETTINGS);
        intent.setData(Uri.parse("package:" + getPackageName()));
        startActivity(intent);
    }

    public static void testme() {
//        double totalTime = 0;
//        Constants.WaitForFeedbackTime = Constants.WaitForFeedbackTimeDefault + Constants.SyncLag;
//        Constants.WaitForSoundingTime = Constants.WaitForSoundingTimeDefault + Constants.SyncLag - Constants.SoundingOffset;
//        Constants.WaitForBerTime = Constants.WaitForBerTimeDefault + Constants.SyncLag;
//        Constants.WaitForPerTime = Constants.WaitForPerTimeDefault + Constants.SyncLag;
//        Constants.WaitForDataTime = Constants.WaitForPerTime;
//
//        double soundingTimeTx = 1;
//        double extractionFeedbackTime = 1;
//        totalTime += (soundingTimeTx + Constants.WaitForFeedbackTime +
//                extractionFeedbackTime);
//        if (Constants.SEND_DATA) {
//            totalTime+=Constants.WaitForDataTime;
//        }
//        Constants.AliceTime = (int)totalTime;
//
//        totalTime=0;
//        int extractSoundingTime = 1;
//        int sendFeedbackTime = 1;
//        totalTime += Constants.WaitForSoundingTime+
//                extractSoundingTime+sendFeedbackTime;
//        totalTime += Constants.SoundingOffset;
//        if (Constants.SEND_DATA) {
//            totalTime+=Constants.WaitForDataTime;
//        }
//        Constants.BobTime = (int)totalTime;



//        totalTime *= 1000;
//        Constants.BobTime = totalTime/1000;
//        //5 symbols, 2 bins
//        double [][][] spect = {
//                {{-1,2,-3,2,1.5}, {1,2,2,1, -1}}, //real
//                {{-2,2,-2,2.5,1.5}, {2,1,1,1,-1 }} //imaginary
//        };
//
//        //5 symbols, 2 bins
//        double [][] gt = {
//                {-1, 1, -1, 1,1}, {1,1,1,1,-1}
//        };
//
//
//        double[] f_seq = FileOperations.readrawasset(av, R.raw.f_seq, 1);
//        double[] SNR = FileOperations.readrawasset(av, R.raw.snr, 1);
//        long t1 = System.currentTimeMillis();
//        double[] test = SNR_freq.calculate_snr(spect, gt);
//        double[] selected = Fre_adaptation.select_fre_bins2(SNR, f_seq, 6);
//        Log.e("time",(System.currentTimeMillis()-t1)+"");
//        select_fre_bins
    }

    public void uiSetup() {
        Constants.gview = (GraphView) findViewById(R.id.graphProd);
        Constants.gview2 = (GraphView) findViewById(R.id.graphProd2);
        Constants.gview3 = (GraphView) findViewById(R.id.graphProd3);
        Constants.sw1 = (Switch) findViewById(R.id.switch1);
        Constants.sw2 = (Switch) findViewById(R.id.switch2);
        Constants.sw3 = (Switch) findViewById(R.id.switch3);
        Constants.sw4 = (Switch) findViewById(R.id.switch4);
        Constants.sw5 = (Switch) findViewById(R.id.switch5);
        Constants.sw6 = (Switch) findViewById(R.id.switch6);
        Constants.sw7 = (Switch) findViewById(R.id.switch7);
        Constants.sw8 = (Switch) findViewById(R.id.switch8);
        Constants.sw9 = (Switch) findViewById(R.id.switch9);
        Constants.sw10 = (Switch) findViewById(R.id.switch10);
        Constants.sw11 = (Switch) findViewById(R.id.switch11);
        Constants.sw12 = (Switch) findViewById(R.id.switch12);
        Constants.sw13 = (Switch) findViewById(R.id.switch13);
        Constants.spinner = (Spinner) findViewById(R.id.spinner);
        Constants.spinner2 = (Spinner) findViewById(R.id.spinner2);
        Constants.spinner3 = (Spinner) findViewById(R.id.spinner3);

        ArrayList<String> arrayList = new ArrayList<>();
        arrayList.add("1/2");
        arrayList.add("2/3");

        ArrayAdapter<String> arrayAdapter = new ArrayAdapter<String>(this,
                android.R.layout.simple_spinner_item, arrayList);
        arrayAdapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        Constants.spinner.setAdapter(arrayAdapter);
        Constants.spinner.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {
            @Override
            public void onItemSelected(AdapterView<?> parent, View view, int position, long id) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                if (position == 0) {
                    editor.putString("code_rate", Constants.CodeRate.C1_2.toString());
                    Constants.codeRate = Constants.CodeRate.C1_2;
                }
                else if (position == 1) {
                    editor.putString("code_rate", Constants.CodeRate.C2_3.toString());
                    Constants.codeRate = Constants.CodeRate.C2_3;
                }
                editor.commit();
            }
            @Override
            public void onNothingSelected(AdapterView <?> parent) {
            }
        });

        ArrayList<String> arrayList2 = new ArrayList<>();
        arrayList2.add("1");
        arrayList2.add("2");
        ArrayAdapter<String> arrayAdapter2 = new ArrayAdapter<String>(this,
                android.R.layout.simple_spinner_item, arrayList2);
        arrayAdapter2.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        Constants.spinner2.setAdapter(arrayAdapter2);
        Constants.spinner2.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {
            @Override
            public void onItemSelected(AdapterView<?> parent, View view, int position, long id) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                Log.e("snr","pos "+position+"");
                if (position == 0) {
                    editor.putInt("snr_method", 1);
                    Constants.snr_method = 1;
                    Constants.et4.setText(Constants.SNR_THRESH1+"");
                }
                else if (position == 1) {
                    editor.putInt("snr_method", 2);
                    Constants.snr_method = 2;
                    Constants.et4.setText(Constants.SNR_THRESH2+"");
                }
                editor.commit();
            }
            @Override
            public void onNothingSelected(AdapterView <?> parent) {
            }
        });

        ArrayList<String> arrayList3 = new ArrayList<>();
        arrayList3.add("960");
        arrayList3.add("1920");
        arrayList3.add("4800");
        arrayList3.add("9600");

        ArrayAdapter<String> arrayAdapter3 = new ArrayAdapter<String>(this,
                android.R.layout.simple_spinner_item, arrayList3);
        arrayAdapter3.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        Constants.spinner3.setAdapter(arrayAdapter3);
        Constants.spinner3.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {
            @Override
            public void onItemSelected(AdapterView<?> parent, View view, int position, long id) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                editor.putInt("ns", Integer.parseInt(arrayList3.get(position)));
                Constants.Ns = Integer.parseInt(arrayList3.get(position));
                editor.commit();
                Constants.updateNbins();
            }
            @Override
            public void onNothingSelected(AdapterView <?> parent) {
            }
        });

        Constants.et1 = (EditText) findViewById(R.id.editTextNumber);
        Constants.et2 = (EditText) findViewById(R.id.editTextNumber2);
        Constants.et3 = (EditText) findViewById(R.id.editTextNumber3);
        Constants.et4 = (EditText) findViewById(R.id.editTextNumber4);
        Constants.et5 = (EditText) findViewById(R.id.editTextNumber5);
        Constants.et6 = (EditText) findViewById(R.id.editTextNumber6);
        Constants.et7 = (EditText) findViewById(R.id.editTextNumber7);
        Constants.et8 = (EditText) findViewById(R.id.editTextNumber8);
        Constants.et9 = (EditText) findViewById(R.id.editTextNumber9);
        Constants.et10 = (EditText) findViewById(R.id.editTextNumber10);
        Constants.et11 = (EditText) findViewById(R.id.editTextNumber11);
        Constants.et12 = (EditText) findViewById(R.id.editTextNumber12);
        Constants.et13 = (EditText) findViewById(R.id.editTextNumber13);
        Constants.et14 = (EditText) findViewById(R.id.editTextNumber14);
        Constants.et15 = (EditText) findViewById(R.id.editTextNumber15);
        Constants.et17 = (EditText) findViewById(R.id.editTextNumber17);
        Constants.et18 = (EditText) findViewById(R.id.editTextNumber18);
        Constants.et25 = (EditText) findViewById(R.id.editTextNumber25);
        Constants.et26 = (EditText) findViewById(R.id.editTextNumber26);
        Constants.et27 = (EditText) findViewById(R.id.editTextNumber27);

        Constants.tv1 = (TextView) findViewById(R.id.textView);
        Constants.tv2 = (TextView) findViewById(R.id.textView2);
        Constants.tv3 = (TextView) findViewById(R.id.textView3);
        Constants.tv4 = (TextView) findViewById(R.id.textView4);
        Constants.tv5 = (TextView) findViewById(R.id.textView5);
        Constants.debugPane = (TextView) findViewById(R.id.debugPane);
        Constants.tv7 = (TextView) findViewById(R.id.textView7);
        Constants.tv8 = (TextView) findViewById(R.id.textView8);
        Constants.tv9 = (TextView) findViewById(R.id.textView9);
        Constants.tv10 = (TextView) findViewById(R.id.textView10);
        Constants.tv13 = (TextView) findViewById(R.id.textView13);
        Constants.tv14 = (TextView) findViewById(R.id.textView14);
        Constants.tv15 = (TextView) findViewById(R.id.textView15);
        Constants.tv16 = (TextView) findViewById(R.id.textView16);
        Constants.tv17 = (TextView) findViewById(R.id.textView17);
        Constants.tv18 = (TextView) findViewById(R.id.textView18);
        Constants.tv19 = (TextView) findViewById(R.id.textView19);
        Constants.tv20 = (TextView) findViewById(R.id.textView20);
        Constants.tv21 = (TextView) findViewById(R.id.textView21);

        Constants.sview = (NestedScrollView) findViewById(R.id.scrollView);
        Constants.tv5.setMovementMethod(new ScrollingMovementMethod());
        Constants.startButton = (Button) findViewById(R.id.button);
        Constants.clearButton = (Button) findViewById(R.id.button2);
        Constants.stopButton = (Button) findViewById(R.id.button3);
        Constants.tv6 = (TextView) findViewById(R.id.textView6);
        av = this;

        Constants.sw2.setVisibility(View.GONE);
        Constants.sw3.setVisibility(View.GONE);
        Constants.sw4.setVisibility(View.GONE);
        Constants.sw5.setVisibility(View.GONE);
        Constants.sw6.setVisibility(View.GONE);
        Constants.sw7.setVisibility(View.GONE);
        Constants.sw8.setVisibility(View.GONE);
//        Constants.sw9.setVisibility(View.GONE);

        Constants.et2.setVisibility(View.GONE);
//        Constants.et4.setVisibility(View.GONE);
        Constants.et5.setVisibility(View.GONE);
        Constants.et6.setVisibility(View.GONE);
        Constants.et7.setVisibility(View.GONE);
        Constants.et8.setVisibility(View.GONE);
        Constants.et11.setVisibility(View.GONE);
//        Constants.et13.setVisibility(View.GONE);
        Constants.et14.setVisibility(View.GONE);
        Constants.et15.setVisibility(View.GONE);
        Constants.et17.setVisibility(View.GONE);
        Constants.et18.setVisibility(View.GONE);

        Constants.tv2.setVisibility(View.GONE);
//        Constants.tv5.setVisibility(View.GONE);
        Constants.tv7.setVisibility(View.GONE);
        Constants.tv8.setVisibility(View.GONE);
        Constants.tv9.setVisibility(View.GONE);
        Constants.tv10.setVisibility(View.GONE);
        Constants.tv13.setVisibility(View.GONE);
//        Constants.tv14.setVisibility(View.GONE);
        Constants.tv15.setVisibility(View.GONE);
//        Constants.tv17.setVisibility(View.GONE);
        Constants.tv18.setVisibility(View.GONE);
        Constants.tv19.setVisibility(View.GONE);
        Constants.tv20.setVisibility(View.GONE);
        Constants.tv21.setVisibility(View.GONE);

        Constants.spinner2.setVisibility(View.GONE);

        Constants.sw1.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                if (isChecked) {
                    editor.putString("user", Constants.User.Alice.toString());
                    Constants.user  = Constants.User.Alice;
                    Constants.sw1.setText(Constants.User.Alice.toString());
                    Constants.sw2.setEnabled(false);
                }
                else {
                    editor.putString("user", Constants.User.Bob.toString());
                    Constants.user  = Constants.User.Bob;
                    Constants.sw1.setText(Constants.User.Bob.toString());
                    Constants.sw2.setEnabled(true);
                }
                editor.commit();
            }
        });

        Constants.sw2.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                editor.putBoolean("decode_data", isChecked);
                Constants.DecodeData  = isChecked;
                editor.commit();
            }
        });

        Constants.sw3.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                editor.putBoolean("test", isChecked);
                Constants.TEST  = isChecked;
                editor.commit();
            }
        });

        Constants.sw4.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                editor.putBoolean("imu", isChecked);
                Constants.imu  = isChecked;
                editor.commit();
            }
        });

        Constants.sw5.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                editor.putBoolean("stereo", isChecked);
                Constants.stereo  = isChecked;
                editor.commit();
            }
        });

        Constants.sw6.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                if (isChecked) {
                    editor.putString("est_sig", Constants.EstSignalType.Chirp.toString());
                    Constants.est_sig = Constants.EstSignalType.Chirp;
                }
                else {
                    editor.putString("est_sig", Constants.EstSignalType.Symbol.toString());
                    Constants.est_sig = Constants.EstSignalType.Symbol;
                }
                Constants.sw6.setText(Constants.est_sig.toString());
                editor.commit();
            }
        });

        Constants.sw7.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                editor.putBoolean("feed_pre", isChecked);
                Constants.feedbackPreamble  = isChecked;
                editor.commit();
            }
        });

        Constants.sw8.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                editor.putBoolean("send_data", isChecked);
                Constants.SEND_DATA  = isChecked;
                editor.commit();
            }
        });

        Constants.sw9.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                editor.putBoolean("flip_symbol", isChecked);
                Constants.FLIP_SYMBOL  = isChecked;
                if (isChecked) {
                    Constants.normal_syms = new int[]{1,3,5,7};
                }
                else {
                    Constants.normal_syms = new int[]{0,1,2,3,4,5,6,7};
                }
                editor.commit();
            }
        });

        Constants.sw10.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                editor.putBoolean("naiser", isChecked);
                Constants.NAISER  = isChecked;
                Constants.updateNaiser(MainActivity.av);
                editor.commit();
            }
        });

        Constants.sw11.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                editor.putBoolean("check_sym", isChecked);
                Constants.CHECK_SYM  = isChecked;
                editor.commit();
            }
        });

        Constants.sw12.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                editor.putBoolean("data_naiser", isChecked);
                Constants.DATA_NAISER  = isChecked;
                editor.commit();
            }
        });

        Constants.sw13.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                editor.putBoolean("onebin", isChecked);
                Constants.onebin  = isChecked;
                editor.commit();
            }
        });

        Constants.et1.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et1.getText().toString();
                if (Utils.isFloat(ss)) {
                    try {
                        editor.putFloat("volume", Float.parseFloat(ss));
                        editor.commit();
                        Constants.volume = Float.parseFloat(ss);
                    }
                    catch(Exception e) {
                        Log.e(Constants.LOG,e.getMessage());
                    }
                }
            }
        });

//        Constants.et2.addTextChangedListener(new TextWatcher() {
//            @Override
//            public void afterTextChanged(Editable s) {
//            }
//            @Override
//            public void beforeTextChanged(CharSequence s, int start,
//                                          int count, int after) {
//            }
//            @Override
//            public void onTextChanged(CharSequence s, int start,
//                                      int before, int count) {
//                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
//                String ss = Constants.et2.getText().toString();
//                if (Utils.isInteger(ss)) {
//                    editor.putInt("preamble_len", Integer.parseInt(ss));
//                    editor.commit();
//                    Constants.preambleTime = Integer.parseInt(ss);
//                }
//            }
//        });

        Constants.et3.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et3.getText().toString();
                if (Utils.isInteger(ss)) {
                    editor.putInt("init_sleep", Integer.parseInt(ss));
                    editor.commit();
                    Constants.initSleep = Integer.parseInt(ss);
                }
            }
        });

        Constants.et4.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et4.getText().toString();
                if (Utils.isInteger(ss)) {
//                    if (Constants.snr_method==1) {
//                        editor.putInt("snr_thresh1", Integer.parseInt(ss));
//                        Constants.SNR_THRESH1 = Integer.parseInt(ss);
//                    }
//                    else if (Constants.snr_method==2) {
                        editor.putInt("snr_thresh2", Integer.parseInt(ss));
                        Constants.SNR_THRESH2 = Integer.parseInt(ss);
//                    }
                    editor.commit();
                }
            }
        });

        Constants.et5.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et5.getText().toString();
                if (Utils.isInteger(ss)) {
                    editor.putInt("nsyms", Integer.parseInt(ss));
                    editor.commit();
                    Constants.Nsyms = Integer.parseInt(ss);
                }
            }
        });

        Constants.et6.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et6.getText().toString();
                if (Utils.isInteger(ss)) {
                    editor.putInt("f1", Integer.parseInt(ss));
                    editor.commit();
                    Constants.f_range[0] = Integer.parseInt(ss);
                    Constants.updateNbins();
                }
            }
        });

        Constants.et7.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et7.getText().toString();
                if (Utils.isInteger(ss)) {
                    editor.putInt("f2", Integer.parseInt(ss));
                    editor.commit();
                    Constants.f_range[1] = Integer.parseInt(ss);
                    Constants.updateNbins();
                }
            }
        });

        Constants.et8.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et8.getText().toString();
                if (Utils.isInteger(ss)) {
                    editor.putInt("symreps", Integer.parseInt(ss));
                    editor.commit();
                    Constants.data_symreps = Integer.parseInt(ss);
                }
            }
        });

        Constants.et9.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et9.getText().toString();
                if (Utils.isInteger(ss)) {
                    editor.putInt("mattempts", Integer.parseInt(ss));
                    editor.commit();
                    Constants.mattempts = Integer.parseInt(ss);
                }
            }
        });

        Constants.et10.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et10.getText().toString();
                if (Utils.isInteger(ss)) {
                    editor.putInt("exp_num", Integer.parseInt(ss));
                    editor.commit();
                    Constants.exp_num = Integer.parseInt(ss);
                }
            }
        });
        Constants.et12.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et12.getText().toString();
                if (Utils.isFloat(ss)) {
                    editor.putFloat("scale_factor", Float.parseFloat(ss));
                    editor.commit();
                    Constants.FreAdaptScaleFactor = Float.parseFloat(ss);
                }
            }
        });
        Constants.et13.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et13.getText().toString();
                if (Utils.isInteger(ss)) {
                    editor.putInt("snr_thresh2_2", Integer.parseInt(ss));
                    editor.commit();
                    Constants.SNR_THRESH2_2 = Integer.parseInt(ss);
                }
            }
        });

        Constants.et14.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et14.getText().toString();
                if (Utils.isFloat(ss)) {
                    editor.putFloat("xcorr_thresh", (float)Double.parseDouble(ss));
                    editor.commit();
                    Constants.MinXcorrVal = Double.parseDouble(ss);
                }
            }
        });

        Constants.et15.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et15.getText().toString();
                if (Utils.isFloat(ss)) {
                    editor.putFloat("xcorr_thresh2", (float)Double.parseDouble(ss));
                    editor.commit();
                    Constants.XCORR_MAX_VAL_HEIGHT_FAC = Double.parseDouble(ss);
                }
            }
        });

        Constants.et17.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et17.getText().toString();
                if (Utils.isInteger(ss)) {
                    editor.putInt("var_thresh", Integer.parseInt(ss));
                    editor.commit();
                    Constants.VAR_THRESH = Integer.parseInt(ss);
                }
            }
        });
        Constants.et18.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et18.getText().toString();
                if (Utils.isInteger(ss)) {
                    editor.putInt("xcorr_above_thresh", Integer.parseInt(ss));
                    editor.commit();
                    Constants.XcorrAboveThresh = Integer.parseInt(ss);
                }
            }
        });

        Constants.et25.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et25.getText().toString();
                if (Utils.isFloat(ss)) {
                    editor.putFloat("naiser_thresh", Float.parseFloat(ss));
                    editor.commit();
                    Constants.NaiserThresh = Float.parseFloat(ss);
                }
            }
        });

        Constants.et26.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et26.getText().toString();
                if (Utils.isInteger(ss)) {
                    editor.putInt("feedback_thresh", Integer.parseInt(ss));
                    editor.commit();
                    Constants.FEEDBACK_SNR_THRESH = Integer.parseInt(ss);
                }
            }
        });

        Constants.et27.addTextChangedListener(new TextWatcher() {
            @Override
            public void afterTextChanged(Editable s) {
            }
            @Override
            public void beforeTextChanged(CharSequence s, int start,
                                          int count, int after) {
            }
            @Override
            public void onTextChanged(CharSequence s, int start,
                                      int before, int count) {
                SharedPreferences.Editor editor = PreferenceManager.getDefaultSharedPreferences(av).edit();
                String ss = Constants.et27.getText().toString();
                if (Utils.isInteger(ss)) {
                    editor.putInt("checksym_snrthresh", Integer.parseInt(ss));
                    editor.commit();
                    Constants.CheckSymSNRThresh = Integer.parseInt(ss);
                }
            }
        });
    }

    protected void onPause() {
        super.onPause();
        sensorManager.unregisterListener(this);
    }

    protected void onResume() {
        super.onResume();
        sensorManager.registerListener(this, accelerometer, SensorManager.SENSOR_DELAY_FASTEST);
        sensorManager.registerListener(this, gyroscope, SensorManager.SENSOR_DELAY_FASTEST);
        FullScreencall();
    }

    @Override
    public void onSensorChanged(SensorEvent event) {
        if (Constants.sensorFlag && Constants.imu) {
            if (event.sensor.equals(accelerometer)) {
                Constants.acc.add(event.values[0]+","+event.values[1]+","+event.values[2]+"\n");
            }
            else if (event.sensor.equals(gyroscope)) {
                Constants.gyro.add(event.values[0]+","+event.values[1]+","+event.values[2]+"\n");
            }
        }
    }

    @Override
    public void onAccuracyChanged(Sensor sensor, int i) {

    }

    public void onstart(View v) {
        sensorManager.registerListener(this, accelerometer, SensorManager.SENSOR_DELAY_FASTEST);
        sensorManager.registerListener(this, gyroscope, SensorManager.SENSOR_DELAY_FASTEST);

        Constants.ts = System.currentTimeMillis();

        Constants.work=true;
        FileOperations.mkdir(av,Utils.getDirName());
//        FileOperations.writetofile(av, Constants.ts+"", Utils.genName(Constants.SignalType.Timestamp,0)+".txt");

        Constants.tv6.setText(Utils.trimmed_ts());
        Constants.task = new SendChirpAsyncTask(this,Constants.mattempts);
        Constants.task.executeOnExecutor(AsyncTask.THREAD_POOL_EXECUTOR);
    }

    public void onstop(View v) {
        if (Constants.task != null) {
            Constants.task.cancel(true);
        }
        if (Constants.timer!=null) {
            Constants.timer.cancel();
        }

        sensorManager.unregisterListener(this);
        Constants.work=false;
        Log.e("asdf","onstop");
        Constants.sensorFlag=false;
        if (Constants.acc != null && Constants.acc.size() > 0) {
            FileOperations.writeSensors(this, Constants.ts+".txt");
        }
        if (Constants._OfflineRecorder!=null) {
            Constants._OfflineRecorder.halt2();
        }
        if (Constants.sp1!=null && Constants.sp1.track1!=null&&
            Constants.sp1.track1.getState()== AudioTrack.STATE_INITIALIZED) {
            Constants.sp1.pause();
        }
        Constants.toggleUI(true);
    }

    @Override
    public void onWindowFocusChanged(boolean hasFocus) {
        super.onWindowFocusChanged(hasFocus);
        FullScreencall();
    }

    public void FullScreencall() {
        if(Build.VERSION.SDK_INT > 11 && Build.VERSION.SDK_INT < 19) { // lower api
            View v = this.getWindow().getDecorView();
            v.setSystemUiVisibility(View.GONE);
        } else if(Build.VERSION.SDK_INT >= 19) {
            //for new api versions.
            View decorView = getWindow().getDecorView();
            int uiOptions = View.SYSTEM_UI_FLAG_HIDE_NAVIGATION | View.SYSTEM_UI_FLAG_IMMERSIVE_STICKY;
            decorView.setSystemUiVisibility(uiOptions);
        }
    }

    public void clear(View v) {
        Constants.debugPane.setText("");
        FullScreencall();
    }

    public static void unreg(Activity av) {
        sensorManager.unregisterListener((MainActivity)av);
    }

}
