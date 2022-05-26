#include <jni.h>
#include <string.h>
#include <fftw3.h>
#include <math.h>

extern "C"
JNIEXPORT jdoubleArray JNICALL
Java_com_example_root_ffttest3_Utils_fftnative_1double(JNIEnv *env, jobject thiz, jdoubleArray data,
                                                         jint N) {
    fftw_complex *in , *out;
    fftw_plan p;

    jdouble *doubleArray = env->GetDoubleArrayElements(data, NULL);
    int datalen = env -> GetArrayLength(data);

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * datalen);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * datalen);

    for (int i = 0; i < datalen; i++) {
        in[i][0] = 0;
        in[i][1] = 0;
        out[i][0] = 0;
        out[i][1] = 0;
    }

    for (int i = 0; i < datalen; i++) {
        in[i][0] = doubleArray[i];
    }

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    jdouble mag[N];
    for (int i = 0; i < N; i++) {
        double real = out[i][0];
        double imag = out[i][1];

        mag[i] = sqrt((real*real)+(imag*imag));
//        mag[i] = 20*log10(mag[i]);
    }

    jdoubleArray result;
    result = env->NewDoubleArray(N);
    env->SetDoubleArrayRegion(result, 0, N, mag);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    return result;
}

extern "C"
JNIEXPORT jdoubleArray JNICALL
Java_com_example_root_ffttest3_Utils_fftnative_1short(JNIEnv *env, jobject thiz, jshortArray data,
                                                        jint N) {
    fftw_complex *in , *out;
    fftw_plan p;

    jshort *shortArray = env->GetShortArrayElements(data, NULL);
    int datalen = env -> GetArrayLength(data);

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * datalen);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * datalen);

    for (int i = 0; i < datalen; i++) {
        in[i][0] = 0;
        in[i][1] = 0;
        out[i][0] = 0;
        out[i][1] = 0;
    }

    for (int i = 0; i < datalen; i++) {
        in[i][0] = shortArray[i];
    }

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    jdouble mag[N];
    for (int i = 0; i < N; i++) {
        double real = out[i][0];
        double imag = out[i][1];

        mag[i] = sqrt((real*real)+(imag*imag));
//        mag[i] = 20*log10(mag[i]);
    }

    jdoubleArray result;
    result = env->NewDoubleArray(N);
    env->SetDoubleArrayRegion(result, 0, N, mag);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    return result;
}

extern "C"
JNIEXPORT jobjectArray JNICALL
Java_com_example_root_ffttest3_Utils_fftcomplexinoutnative_1double(JNIEnv *env, jobject thiz,
                                                                 jobjectArray data, jint N) {

    jdoubleArray realar1 = (jdoubleArray) env->GetObjectArrayElement(data, 0);
    jdoubleArray imagar1 = (jdoubleArray) env->GetObjectArrayElement(data, 1);

    jdouble *real1 = env->GetDoubleArrayElements(realar1, NULL);
    jdouble *imag1 = env->GetDoubleArrayElements(imagar1, NULL);

    jint datalen = env -> GetArrayLength(realar1);

    jint fft_len = N;
    if (datalen > N) {
        fft_len = datalen;
    }

    fftw_complex *in , *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft_len);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft_len);

    for (int i = 0; i < fft_len; i++) {
        in[i][0] = real1[i];
        in[i][1] = imag1[i];
        out[i][0] = 0;
        out[i][1] = 0;
    }

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    jdouble real[N];
    jdouble imag[N];
    for (int i = 0; i < N; i++) {
        real[i] = out[i][0];
        imag[i] = out[i][1];
    }

    jdoubleArray realResult;
    jdoubleArray imagResult;
    realResult = env->NewDoubleArray(N);
    imagResult = env->NewDoubleArray(N);
    env->SetDoubleArrayRegion(realResult, 0, N, real);
    env->SetDoubleArrayRegion(imagResult, 0, N, imag);

    jobjectArray outarray = env->NewObjectArray(2, env->GetObjectClass(realResult), 0);
    env->SetObjectArrayElement(outarray, 0, realResult);
    env->SetObjectArrayElement(outarray, 1, imagResult);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    return outarray;
}

extern "C"
JNIEXPORT jobjectArray JNICALL
Java_com_example_root_ffttest3_Utils_fftcomplexoutnative_1double(JNIEnv *env, jobject thiz,
                                                                   jdoubleArray data, jint N) {

    fftw_complex *in , *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++) {
        in[i][0] = 0;
        in[i][1] = 0;
        out[i][0] = 0;
        out[i][1] = 0;
    }

    jdouble *doubleArray = env->GetDoubleArrayElements(data, NULL);
    int datalen = env -> GetArrayLength(data);
    for (int i = 0; i < datalen; i++) {
        in[i][0] = doubleArray[i];
    }

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    jdouble real[N];
    jdouble imag[N];
    for (int i = 0; i < N; i++) {
        real[i] = out[i][0];
        imag[i] = out[i][1];
    }

    jdoubleArray realResult;
    jdoubleArray imagResult;
    realResult = env->NewDoubleArray(N);
    imagResult = env->NewDoubleArray(N);
    env->SetDoubleArrayRegion(realResult, 0, N, real);
    env->SetDoubleArrayRegion(imagResult, 0, N, imag);

    jobjectArray outarray = env->NewObjectArray(2, env->GetObjectClass(realResult), 0);
    env->SetObjectArrayElement(outarray, 0, realResult);
    env->SetObjectArrayElement(outarray, 1, imagResult);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    return outarray;
}

extern "C"
JNIEXPORT jobjectArray JNICALL
Java_com_example_root_ffttest3_Utils_fftcomplexoutnative_1short(JNIEnv *env, jobject thiz,
                                                                  jshortArray data, jint N) {

    fftw_complex *in , *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++) {
        in[i][0] = 0;
        in[i][1] = 0;
        out[i][0] = 0;
        out[i][1] = 0;
    }

    jshort *shortArray = env->GetShortArrayElements(data, NULL);
    int datalen = env -> GetArrayLength(data);
    for (int i = 0; i < datalen; i++) {
        in[i][0] = shortArray[i];
    }

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    jdouble real[N];
    jdouble imag[N];
    for (int i = 0; i < N; i++) {
        real[i] = out[i][0];
        imag[i] = out[i][1];
    }

    jdoubleArray realResult;
    jdoubleArray imagResult;
    realResult = env->NewDoubleArray(N);
    imagResult = env->NewDoubleArray(N);
    env->SetDoubleArrayRegion(realResult, 0, N, real);
    env->SetDoubleArrayRegion(imagResult, 0, N, imag);

    jobjectArray outarray = env->NewObjectArray(2, env->GetObjectClass(realResult), 0);
    env->SetObjectArrayElement(outarray, 0, realResult);
    env->SetObjectArrayElement(outarray, 1, imagResult);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    return outarray;
}

extern "C"
JNIEXPORT jdoubleArray JNICALL
Java_com_example_root_ffttest3_Utils_ifftnative(JNIEnv *env, jobject thiz,
                                                  jobjectArray data) {
    jdoubleArray real = (jdoubleArray) env->GetObjectArrayElement(data, 0);
    jdoubleArray imag = (jdoubleArray) env->GetObjectArrayElement(data, 1);

    jint N = env -> GetArrayLength(real);

    fftw_complex *in , *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++) {
        in[i][0] = 0;
        in[i][1] = 0;
        out[i][0] = 0;
        out[i][1] = 0;
    }

    jdouble *realArray = env->GetDoubleArrayElements(real, NULL);
    jdouble *imagArray = env->GetDoubleArrayElements(imag, NULL);
    for (int i = 0; i < N; i++) {
        in[i][0] = realArray[i];
    }
    for (int i = 0; i < N; i++) {
        in[i][1] = imagArray[i];
    }

    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    jdouble realout[N/2];
    int counter = 0;
    for (int i = 0; i < N; i+=2) {
        realout[counter++] = out[i][0];
    }

    jdoubleArray result;
    result = env->NewDoubleArray(N/2);
    env->SetDoubleArrayRegion(result, 0, N/2, realout);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    return result;
}

extern "C"
JNIEXPORT jobjectArray JNICALL
Java_com_example_root_ffttest3_Utils_ifftnative2(JNIEnv *env, jobject thiz,
                                                   jobjectArray data) {
    jdoubleArray real = (jdoubleArray) env->GetObjectArrayElement(data, 0);
    jdoubleArray imag = (jdoubleArray) env->GetObjectArrayElement(data, 1);

    jint N = env -> GetArrayLength(real);

    fftw_complex *in , *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++) {
        in[i][0] = 0;
        in[i][1] = 0;
        out[i][0] = 0;
        out[i][1] = 0;
    }

    jdouble *realArray = env->GetDoubleArrayElements(real, NULL);
    jdouble *imagArray = env->GetDoubleArrayElements(imag, NULL);
    for (int i = 0; i < N; i++) {
        in[i][0] = realArray[i];
    }
    for (int i = 0; i < N; i++) {
        in[i][1] = imagArray[i];
    }

    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    jdouble real2[N];
    jdouble imag2[N];
    for (int i = 0; i < N; i++) {
        real2[i] = out[i][0];
        imag2[i] = out[i][1];
    }

    jdoubleArray realResult;
    jdoubleArray imagResult;
    realResult = env->NewDoubleArray(N);
    imagResult = env->NewDoubleArray(N);
    env->SetDoubleArrayRegion(realResult, 0, N, real2);
    env->SetDoubleArrayRegion(imagResult, 0, N, imag2);

    jobjectArray outarray = env->NewObjectArray(2, env->GetObjectClass(realResult), 0);
    env->SetObjectArrayElement(outarray, 0, realResult);
    env->SetObjectArrayElement(outarray, 1, imagResult);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    return outarray;
}

extern "C" JNIEXPORT jobjectArray JNICALL Java_com_example_root_ffttest3_Utils_timesnative(
        JNIEnv *env,
        jclass,
        jobjectArray c1,
        jobjectArray c2) {

    jdoubleArray realar1 = (jdoubleArray) env->GetObjectArrayElement(c1, 0);
    jdoubleArray imagar1 = (jdoubleArray) env->GetObjectArrayElement(c1, 1);
    jdoubleArray realar2 = (jdoubleArray) env->GetObjectArrayElement(c2, 0);
    jdoubleArray imagar2 = (jdoubleArray) env->GetObjectArrayElement(c2, 1);

    jdouble *real1 = env->GetDoubleArrayElements(realar1, NULL);
    jdouble *imag1 = env->GetDoubleArrayElements(imagar1, NULL);
    jdouble *real2 = env->GetDoubleArrayElements(realar2, NULL);
    jdouble *imag2 = env->GetDoubleArrayElements(imagar2, NULL);

    jint N = env -> GetArrayLength(realar1);

    jint counter = 0;
    jdouble real[N];
    jdouble imag[N];
    for (int i = 0; i < N; i++) {
        real[counter] = real1[i]*real2[i]-imag1[i]*imag2[i];
        imag[counter++] = imag1[i]*real2[i]+real1[i]*imag2[i];
    }

    jdoubleArray realResult;
    jdoubleArray imagResult;
    realResult = env->NewDoubleArray(N);
    imagResult = env->NewDoubleArray(N);
    env->SetDoubleArrayRegion(realResult, 0, N, real);
    env->SetDoubleArrayRegion(imagResult, 0, N, imag);

    jobjectArray outarray = env->NewObjectArray(2, env->GetObjectClass(realResult), 0);
    env->SetObjectArrayElement(outarray, 0, realResult);
    env->SetObjectArrayElement(outarray, 1, imagResult);

    return outarray;
}

extern "C" JNIEXPORT jobjectArray JNICALL Java_com_example_root_ffttest3_Utils_dividenative(
        JNIEnv *env,
        jclass,
        jobjectArray c1,
        jobjectArray c2) {

    jdoubleArray realar1 = (jdoubleArray) env->GetObjectArrayElement(c1, 0);
    jdoubleArray imagar1 = (jdoubleArray) env->GetObjectArrayElement(c1, 1);
    jdoubleArray realar2 = (jdoubleArray) env->GetObjectArrayElement(c2, 0);
    jdoubleArray imagar2 = (jdoubleArray) env->GetObjectArrayElement(c2, 1);

    jdouble *real1 = env->GetDoubleArrayElements(realar1, NULL);
    jdouble *imag1 = env->GetDoubleArrayElements(imagar1, NULL);
    jdouble *real2 = env->GetDoubleArrayElements(realar2, NULL);
    jdouble *imag2 = env->GetDoubleArrayElements(imagar2, NULL);

    jint N = env -> GetArrayLength(realar1);

    jint counter = 0;
    jdouble real[N];
    jdouble imag[N];
    for (int i = 0; i < N; i++) {
        jdouble a = real1[i];
        jdouble b = imag1[i];
        jdouble c = real2[i];
        jdouble d = imag2[i];
        real[counter] = (a*c+b*d)/(c*c+d*d);
        imag[counter++] = (b*c-a*d)/(c*c+d*d);
    }

    jdoubleArray realResult;
    jdoubleArray imagResult;
    realResult = env->NewDoubleArray(N);
    imagResult = env->NewDoubleArray(N);
    env->SetDoubleArrayRegion(realResult, 0, N, real);
    env->SetDoubleArrayRegion(imagResult, 0, N, imag);

    jobjectArray outarray = env->NewObjectArray(2, env->GetObjectClass(realResult), 0);
    env->SetObjectArrayElement(outarray, 0, realResult);
    env->SetObjectArrayElement(outarray, 1, imagResult);

    return outarray;
}

extern "C" JNIEXPORT void JNICALL Java_com_example_root_ffttest3_Utils_conjnative(
        JNIEnv *env,
        jclass,
        jobjectArray data) {

    jdoubleArray imag = (jdoubleArray) env->GetObjectArrayElement(data, 1);
    jdouble *imagArray = env->GetDoubleArrayElements(imag, NULL);
    jint N = env -> GetArrayLength(imag);

    for (int i = 0; i < N; i++) {
        imagArray[i] = -imagArray[i];
    }
}

extern "C" JNIEXPORT jdoubleArray JNICALL Java_com_example_root_ffttest3_Utils_fir(
        JNIEnv *env,
        jclass,
        jdoubleArray data,
        jdoubleArray h) {

    jdouble *jdata = env->GetDoubleArrayElements(data, NULL);
    jdouble *jh = env->GetDoubleArrayElements(h, NULL);

    jint lenData = env->GetArrayLength(data);
    jint lenH = env->GetArrayLength(h);

    jint nconv = lenH+lenData-1;

    jdoubleArray out;
    out = env->NewDoubleArray(nconv);

    jdouble temp[nconv];

    for (int i=0; i<nconv; i++) {
        temp[i]=0;
    }
//    for (int i=0; i<nconv; i++) {
//        jint x_start = 0;
//        if (i-lenH+1 > 0) {
//            x_start = i-lenH+1;
//        }
//        jint x_end = i+1;
//        if (lenData<i+1) {
//            x_end=lenData;
//        }
//
//        jint h_start=i;
//        if (lenH-1<i) {
//            h_start = lenH-1;
//        }
//
//        for(int j=x_start; j<x_end; j++) {
//            temp[j] += jh[h_start--] * jdata[j];
//        }
//    }
    for (int n = 0; n < nconv; n++){
        jint kmin, kmax;

        temp[n] = 0;

        kmin = (n >= lenH - 1) ? n - (lenH - 1) : 0;
        kmax = (n < lenData - 1) ? n : lenData - 1;

        for (int k = kmin; k <= kmax; k++) {
            temp[n] += jdata[k] * jh[n - k];
        }
    }

    env->SetDoubleArrayRegion(out, 0, nconv, temp);

    return out;
}

extern "C" JNIEXPORT jdoubleArray JNICALL Java_com_example_root_ffttest3_Utils_bandpass(
        JNIEnv *env,
        jclass,
        jdoubleArray data) {

    jint N = env->GetArrayLength(data);
    jdouble *input = env->GetDoubleArrayElements(data, NULL);

    jint m_numBiquads = 6;

    jdouble a0[] = {1.0,1.0,1.0,1.0,1.0,1.0};
    //butterworth
//    jdouble a1[] = {-1.7676841922606346,-1.9120840971742479,-1.704950180189151,-1.8382672336028814,-1.7092159623343837,-1.7643149301109529};
//    jdouble a2[] = {0.9372642800019002,0.964368159585712,0.8446814267964798,0.8968433783987471,0.8133733502161703,0.8395584361646247};
//    jdouble b0[] = {6.241911798109899E-7,1.0,1.0,1.0,1.0,1.0};
//    jdouble b1[] = {1.2483823596219797E-6,-2.0,2.0,-2.0,2.0,-2.0};
//    jdouble b2[] = {6.241911798109899E-7,1.0,1.0,1.0,1.0,1.0};

//bessel
    jdouble a1[] = {-0.8483132820723579,-1.8379673789503772,-0.9104634235790083,-1.7631232694140826,-0.7921543000958807,-1.9173321405089663};
    jdouble a2[] = {0.2834093153669302,0.8491194750058104,0.21989677363118537,0.7781746366113,0.4696113329720846,0.9252519876229148};
    jdouble b0[] = {0.0013832725795677246,1.0,1.0,1.0,1.0,1.0};
    jdouble b1[] = {0.002766545159135449,-2.0,2.0,-2.0,2.0,-2.0};
    jdouble b2[] = {0.0013832725795677246,1.0,1.0,1.0,1.0,1.0};

    jdouble m_v1[] = {0,0,0,0,0,0};
    jdouble m_v2[] = {0,0,0,0,0,0};

    for (int k = 0; k < N; k++) {
        jdouble in = input[k];

        for (int i = 0; i < m_numBiquads; i++) {

            jdouble w = in - a1[i] * m_v1[i] - a2[i] * m_v2[i];
            in = b0[i] * w + b1[i] * m_v1[i] + b2[i] * m_v2[i];

            m_v2[i] = m_v1[i];
            m_v1[i] = w;

        }

        input[k] = in;
    }

    jdoubleArray out;
    out = env->NewDoubleArray(N);
    env->SetDoubleArrayRegion(out, 0, N, input);

    return out;
}