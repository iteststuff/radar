#include <portaudio.h>
#include <stdio.h>
#include <stdlib.h>
#include "GL/glew.h"
#include "GL/glfw.h"
#include <math.h>
#include <fftw3.h>
//#include <glm/glm.hpp>
//#include <GL/gl.h>
//#include <GL/glu.h> 


#define SAMPLE_RATE 44100
#define FRAMES_PER_BUFFER 512
#define SAMPLE_TIME 0.06
#define CHANNELS 2


typedef struct
{
    int          frameIndex;  /* Index into sample array. */
    int          maxFrameIndex;
    float        *recordedSamples;
}
paTestData;

void InitAudioStreaming(paTestData* data, PaStreamParameters* inputParameters,
			PaStreamParameters* outputParameters);

int recordCallback(const void *input,
		     void *output,
		     unsigned long frameCount,
		     const PaStreamCallbackTimeInfo* timeInfo,
		     PaStreamCallbackFlags statusFlags,
		     void *userData);

int playCallback(const void *input,
		     void *output,
		     unsigned long frameCount,
		     const PaStreamCallbackTimeInfo* timeInfo,
		     PaStreamCallbackFlags statusFlags,
		     void *userData);

void ErrExit(char *errString);
char Prompt();
double DecRound(double x, int d);
void Record(paTestData* data, PaStream* stream, PaStreamParameters* inputParameters);
void Play(paTestData* data, PaStream* stream, PaStreamParameters* outputParameters);
void Doppler(paTestData* data, PaStream* stream, PaStreamParameters* inputParameters);
void DopplerSample(paTestData* data, PaStream* stream,
		   fftw_complex* fft_buff, fftw_plan* plan);
void Range(paTestData* data, PaStream* stream, PaStreamParameters* inputParameters);
void GetHalfPeriod(float* dual_chan_buff, int buff_size, float threshold,
		   int* start, int* stop, int* rising);
void RangeSample(paTestData* data, PaStream* stream,
		 fftw_complex* fft_buff, fftw_plan* plan);
void CreateScope();
void Artist();


void InitAudioStreaming(paTestData* data, PaStreamParameters* inputParameters,
			PaStreamParameters* outputParameters){
  int     totalFrames;
  int     numSamples;
  int     i;
  PaError err = paNoError;


  totalFrames = SAMPLE_RATE * SAMPLE_TIME;
  data->maxFrameIndex = totalFrames;
  data->frameIndex = 0;
  numSamples = totalFrames * CHANNELS;
  data->recordedSamples = (float*) malloc(numSamples * sizeof(float));
  
  if(data->recordedSamples == NULL)
    ErrExit("NOOOOO Memory for You!!!!");

  for(i = 0; i < numSamples; i++)
    data->recordedSamples[i] = 0;

  err = Pa_Initialize();
  if(err != paNoError)
    ErrExit("Cant initialize sound card.");

  inputParameters->device = Pa_GetDefaultInputDevice();
  if(inputParameters->device == paNoDevice)
    ErrExit("Could not get default input device.");

  inputParameters->channelCount = CHANNELS;
  inputParameters->sampleFormat = paFloat32;
  inputParameters->suggestedLatency =
    Pa_GetDeviceInfo(inputParameters->device)->defaultLowInputLatency;
  inputParameters->hostApiSpecificStreamInfo = NULL;
}



int recordCallback(const void *input,
		     void *output,
		     unsigned long frameCount,
		     const PaStreamCallbackTimeInfo* timeInfo,
		     PaStreamCallbackFlags statusFlags,
		     void *userData )
{
  paTestData *data = (paTestData*)userData;
  const float *rptr = (const float*)input;
  float *wptr = &data->recordedSamples[data->frameIndex * CHANNELS];
  long framesToCalc;
  long i;
  int finished;
  unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;

  (void) output;
  (void) timeInfo;
  (void) statusFlags;

  if(framesLeft < frameCount){
    framesToCalc = framesLeft;
    finished = paComplete;
  }
  else{
    framesToCalc = frameCount;
    finished = paContinue;
  }
  
  if(input == NULL){
    for(i = 0; i < framesToCalc; i++){
      *wptr++ = 0.0f;
      if(CHANNELS == 2) *wptr++ = 0.0f;
    }
  }
  else{
    for(i = 0; i < framesToCalc; i++){
      *wptr++ = *rptr++;
      if(CHANNELS == 2) *wptr++ = *rptr++;
    }
  }
  data->frameIndex += framesToCalc;
  return finished;
}


int playCallback(const void *input,
		     void *output,
		     unsigned long frameCount,
		     const PaStreamCallbackTimeInfo* timeInfo,
		     PaStreamCallbackFlags statusFlags,
		     void *userData )
{
  paTestData *data = (paTestData*)userData;
  float *rptr = &data->recordedSamples[data->frameIndex * CHANNELS];
  float *wptr = (float*)output;
  unsigned int i;
  int finished;
  unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;

  (void) input;
  (void) timeInfo;
  (void) statusFlags;

  if(framesLeft < frameCount){
    for(i = 0; i < framesLeft; i++){
      *wptr++ = *rptr++;
      if(CHANNELS == 2) *wptr++ = *rptr++;
    }
    for(; i < frameCount; i++){
      *wptr++ = 0;
      if(CHANNELS == 2) *wptr++ = 0;
    }
    data->frameIndex += framesLeft;
    finished = paComplete;
  }
  else{
    for(i = 0; i < frameCount; i++){
      *wptr++ = *rptr++;
      if(CHANNELS == 2) *wptr++ = *rptr++;
    }
    data->frameIndex += frameCount;
    finished = paContinue;
  }
  return finished;
}


void ErrExit(char *errString)
{
  Pa_Terminate();
  printf("Something went wrong: %s Exiting!\n", errString);
  exit(0);
}


char Prompt(){
  char c;

  printf("\nWhat would you like to do?\n");
  printf("R. Detect Range\n");
  printf("D. Get some doppler!\n");
  printf("W. Record some audio.\n");
  printf("P. Play some audio.\n");
  printf("S. Speed Calc.\n");
  printf("A. Artist the window.\n");
  printf("Q. Quit.\n");
  
  c = getchar();
  if(c != '\n')
    getchar();
  
  return c;
}


double DecRound(double x, int d){
  int order = pow(10, d);

  x = x * order;
  return round(x)/order;
}


void Record(paTestData* data, PaStream* stream, PaStreamParameters* inputParameters){
  PaError err;

  data->frameIndex = 0;
  err = Pa_OpenStream(&stream, inputParameters, NULL, SAMPLE_RATE,
		      FRAMES_PER_BUFFER, paClipOff, recordCallback,
		      data);
  if(err != paNoError)
    ErrExit("There was an error opening the stream.");
  
  err = Pa_StartStream(stream);
  if(err != paNoError)
    ErrExit("There was an error starting the stream.");
  
  printf("\n=== Now recording!! ===\n");
  
  while((err = Pa_IsStreamActive(stream)) == 1){
    Pa_Sleep(1000);
    printf("index = %d\n", data->frameIndex);
  }
  if(err < 0)
    ErrExit("Active stream died!");
  
  err = Pa_CloseStream(stream);
  if(err != paNoError)
    ErrExit("Could not close stream.");
  
  printf("Done recording.\n\n");
}


void Play(paTestData* data, PaStream* stream, PaStreamParameters* outputParameters){
  PaError err;
  int i;
  FILE *file;

  fftw_complex *fft_buff;
  fftw_plan plan;
  int N = SAMPLE_RATE * SAMPLE_TIME;
  int max_index = 0;
  float max_value = 0;
  float complex_mag = 0;

  fft_buff = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  for(i = 0; i < N; i++){
    fft_buff[i][0] = data->recordedSamples[2*i+1];
    fft_buff[i][1] = 0;
  }

  plan = fftw_plan_dft_1d(N, fft_buff, fft_buff, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  
  file = fopen("FFT_Values.txt", "w");
  fprintf(file, "Frequency Magnitude\n");
  for(i = 0; i <= (N / 2); i++){
    //printf("%d %f %f\n", i, fft_buff[i][0], fft_buff[i][1]);
    //    fprintf(file, "%d %f %f %f\n", i, 
    //	    fft_buff[i][0], fft_buff[i][1], 
    //	    pow(fft_buff[i][0]*fft_buff[i][0] +
    //		fft_buff[i][1]*fft_buff[i][1], 0.5));
    
    complex_mag = sqrt(fft_buff[i][0]*fft_buff[i][0] +
		       fft_buff[i][1]*fft_buff[i][1]);

    fprintf(file, "%f %f\n", ((float)i * SAMPLE_RATE / N), complex_mag);

    if(complex_mag > max_value){
      max_index = i;
      max_value = complex_mag;
    }
  }
  fclose(file);

  printf("The dominant frequency is: %f\n", ((float)max_index * SAMPLE_RATE / N));

  data->frameIndex = 0;

  outputParameters->device = Pa_GetDefaultOutputDevice();
  if(outputParameters->device == paNoDevice){
    ErrExit("Could not get default output device.");
  }
  
  outputParameters->channelCount = CHANNELS;
  outputParameters->sampleFormat = paFloat32;
  outputParameters->suggestedLatency = 
    Pa_GetDeviceInfo(outputParameters->device)->defaultLowOutputLatency;
  outputParameters->hostApiSpecificStreamInfo = NULL;
  
  printf("\n=== Now playing back!! ===\n");
  err = Pa_OpenStream(&stream, NULL, outputParameters, SAMPLE_RATE,
		      FRAMES_PER_BUFFER, paClipOff, playCallback, data);
  if(err != paNoError)
    ErrExit("There was an error opening the stream.");
  
  if(stream){
    err = Pa_StartStream(stream);
    if(err != paNoError)
      ErrExit("There was an error starting the stream.");
    
    printf("Waiting for playback to finish.\n");
    
    while((err = Pa_IsStreamActive(stream)) == 1) Pa_Sleep(100);
    if(err < 0)
      ErrExit("Active stream died!");
    
    err = Pa_CloseStream(stream);
    if(err != paNoError)
      ErrExit("Could not close stream.");
    
    printf("Done playing sound.\n\n");

    fftw_destroy_plan(plan);
    fftw_free(fft_buff);
  }
}


void Doppler(paTestData* data, PaStream* stream, PaStreamParameters* inputParameters){
  PaError err;
  int i;
  fftw_complex *fft_buff;
  fftw_plan plan;
  int N = SAMPLE_RATE * SAMPLE_TIME;

  data->frameIndex = 0;
  err = Pa_OpenStream(&stream, inputParameters, NULL, SAMPLE_RATE,
		      FRAMES_PER_BUFFER, paClipOff, recordCallback,
		      data);
  if(err != paNoError)
    ErrExit("There was an error opening the stream.");
  
  err = Pa_StartStream(stream);
  if(err != paNoError)
    ErrExit("There was an error starting the stream.");
  
  while((err = Pa_IsStreamActive(stream)) == 1){
    Pa_Sleep(SAMPLE_TIME * 100);
  }

  if(err < 0)
    ErrExit("Active stream died!");
  
  fft_buff = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  for(i = 0; i < N; i++){
    fft_buff[i][0] = data->recordedSamples[2*i+1];
    fft_buff[i][1] = 0;
  }

  plan = fftw_plan_dft_1d(N, fft_buff, fft_buff, FFTW_FORWARD, FFTW_ESTIMATE);
  
  for(;;){
    DopplerSample(data, stream, fft_buff, &plan);
  }

  err = Pa_CloseStream(stream);
  if(err != paNoError)
    ErrExit("Could not close stream.");

  fftw_destroy_plan(plan);
  fftw_free(fft_buff);
}

void DopplerSample(paTestData* data, PaStream* stream,
		   fftw_complex* fft_buff, fftw_plan* plan){
  PaError err;
  int i;
  int N = SAMPLE_RATE * SAMPLE_TIME;
  int max_index = 0;
  float max_value = 0;
  float complex_mag = 0;
  int        c = 299792458;
  float      Fd;
  float      Ft = 2427e6; //Vt=2.35V
  float      Vr;
  float      MPH;

  err = Pa_StopStream(stream);
  data->frameIndex = 0;
  err = Pa_StartStream(stream);
  while((err = Pa_IsStreamActive(stream)) == 1){
    Pa_Sleep(SAMPLE_TIME * 100);
  }

  for(i = 0; i < N; i++){
    fft_buff[i][0] = data->recordedSamples[2*i+1];
    fft_buff[i][1] = 0;
  }

  fftw_execute((*plan));
  
  for(i = 0; i <= (N / 2); i++){   
    complex_mag = sqrt(fft_buff[i][0]*fft_buff[i][0] +
		       fft_buff[i][1]*fft_buff[i][1]);

    if((complex_mag > max_value) && (complex_mag > 20)){
      max_index = i;
      max_value = complex_mag;
    }
  }

  Fd = ((float)max_index * SAMPLE_RATE / N);
  Vr = (Fd * c) / (2 * Ft); // meters/sec
  MPH = (Vr * 3600) / (0.0254 * 12 * 5280);

  //  printf("Speed = %.2f MPH. Frequency = %.2f Hz. Amplitude = %.2f\n", 
  //	 DecRound(MPH, 2), DecRound(Fd, 2), DecRound(max_value, 2));
  printf("Speed = %.2f MPH\n", DecRound(MPH, 2));

}


void Range(paTestData* data, PaStream* stream, PaStreamParameters* inputParameters){
  //  int fstart = 2401; //Vt=2.00V
  //  int fstop  = 2496; //Vt=3.40V
  //  int bw     = fstop - fstart; //95MHz
  float tp     = 20e-3; //pulse time period
  int pulse_size = SAMPLE_RATE * tp;
  //  int c      = 299792458;
  // float rr     = c / (2 * bw * 1e6);
  // float max_range = rr * pulse_size / 2;

  PaError err;
  int i;
  fftw_complex *fft_buff;
  fftw_plan plan;

  int fft_bin = SAMPLE_RATE * SAMPLE_TIME / 3;
  int pulse_start = 0;
  int threshold = 0;

  CreateScope();
  data->frameIndex = 0;
  err = Pa_OpenStream(&stream, inputParameters, NULL, SAMPLE_RATE,
		      FRAMES_PER_BUFFER, paClipOff, recordCallback,
		      data);
  if(err != paNoError)
    ErrExit("There was an error opening the stream.");
  
  err = Pa_StartStream(stream);
  if(err != paNoError)
    ErrExit("There was an error starting the stream.");
  
  while((err = Pa_IsStreamActive(stream)) == 1){
    Pa_Sleep(SAMPLE_TIME * 100);
  }

  if(err < 0)
    ErrExit("Active stream died!");

  fft_buff = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pulse_size);
    
  for(i = 1; i < fft_bin; i++){
    if (data->recordedSamples[2*i] > threshold 
	&& data->recordedSamples[2*(i-1)] < threshold){
      pulse_start = i;
      break;
    }
  }
  printf("fft_bin: %d pulse_start+pulse_size: %d\n", fft_bin, pulse_start + pulse_size);
  //if ((fft_bin - (pulse_start + pulse_size)) < 0)
  //  ErrExit("Unable to obtain full pulse period.");

  for(i = 0; i < pulse_size; i++){
      fft_buff[i][0] = data->recordedSamples[pulse_start + 2*i+1];
      fft_buff[i][1] = 0;
  }
  
  plan = fftw_plan_dft_1d(pulse_size, fft_buff, fft_buff, FFTW_FORWARD, FFTW_ESTIMATE);
  
  for(;;){
    RangeSample(data, stream, fft_buff, &plan);
  }

  err = Pa_CloseStream(stream);
  if(err != paNoError)
    ErrExit("Could not close stream.");

  fftw_destroy_plan(plan);
  fftw_free(fft_buff);
}

void GetHalfPeriod(float* dual_chan_buff, int buff_size, float threshold,
		   int* start, int* stop, int* rising)
{
  int i;
  *rising = -1;
  *stop = -1;

  for(i = *start + 2; i < buff_size; i += 2){
    if(dual_chan_buff[i-2] < threshold && 
       dual_chan_buff[i] > threshold){
      *start = i;
      *rising = 1;
      break;
    }
    else if(dual_chan_buff[i-2] > threshold && 
	    dual_chan_buff[i] < threshold){
      *start = i;
      *rising = 0;
      break;
    }
  }

  if(*rising == 1){
    for(i = *start; i < buff_size; i += 2){
      if(dual_chan_buff[i] < threshold
	 && dual_chan_buff[i-2] > threshold){
	*stop = i;
	break;
      }
    }
  }
  else{
    for(i = *start; i < buff_size; i += 2){
      if(dual_chan_buff[i] > threshold
	 && dual_chan_buff[i-2] < threshold){
	*stop = i;
	break;
      }
    }
  }

  if(*rising == -1 || *stop == -1){
    FILE* fout = fopen("fft_dump.csv", "w");
    for(i=0; i < buff_size; i += 2)
      fprintf(fout, "%f\n", dual_chan_buff[i+1]);
    ErrExit("Unable to obtain full pulse period.");
  } 
}

void RangeSample(paTestData* data, PaStream* stream,
		   fftw_complex* fft_buff, fftw_plan* plan){
  int fstart = 2401; //Vt=2.00V
  int fstop  = 2496; //Vt=3.40V
  int bw     = fstop - fstart; //95MHz
  float tp    = 20e-3; //pulse time period
  int pulse_size = tp * SAMPLE_RATE;
  int c      = 299792458;
  float rr     = c / (2 * bw * 1e6);
  float max_range = rr * pulse_size / 2;

  PaError err;
  int   i;
  int   sample_size = SAMPLE_RATE * SAMPLE_TIME * 2;
  int   fft_bin = SAMPLE_RATE * SAMPLE_TIME / 3;

  int   max_index = 0;
  float max_value = 0;
  float complex_mag = 0;
  float Fr;
  float distance;

  int pulse1_start = 0, pulse1_stop = 0, pulse1_rising = 0;
  int pulse2_start = 0, pulse2_stop = 0 , pulse2_rising = 0;
  float pulse1_val, pulse2_val;
  float threshold = 0.0f;
  FILE* fout;

  err = Pa_StopStream(stream);
  data->frameIndex = 0;
  err = Pa_StartStream(stream);
  while((err = Pa_IsStreamActive(stream)) == 1){
    Pa_Sleep(SAMPLE_TIME * 100);
  }

  GetHalfPeriod(data->recordedSamples, sample_size, threshold,
		&pulse1_start, &pulse1_stop, &pulse1_rising);
  pulse2_start = pulse1_stop -2;
  GetHalfPeriod(data->recordedSamples, sample_size, threshold,
		&pulse2_start, &pulse2_stop, &pulse2_rising);
  if(pulse1_rising != pulse2_rising){
    pulse2_start = pulse2_stop - 2;
    GetHalfPeriod(data->recordedSamples, sample_size, threshold,
		&pulse2_start, &pulse2_stop, &pulse2_rising);
    printf("rising %d %d\n", pulse1_rising, pulse2_rising);
  }
  if(pulse1_rising != pulse2_rising)
    printf("Fuckall\n");
  
  //while((fout = fopen("fuckall.csv", "w")) == NULL);
  for(i = 0; i < fft_bin; i++){
    pulse1_val = data->recordedSamples[pulse1_start + i*2 + 1];
    pulse2_val = data->recordedSamples[pulse2_start + i*2 + 1];
    fft_buff[i][0] = pulse1_val - pulse2_val;
    fft_buff[i][1] = 0;
    
    //  fprintf(fout, "%d, %f, %f\n", i, pulse1_val, pulse2_val);
  }
 
  printf("%d %d\n", pulse1_start, pulse2_start);
  fclose(fout);
  fout = fopen("func_gen.csv", "w");
  for(i = 0; i < sample_size/2; i++)
    fprintf(fout, "%d, %f\n", i, data->recordedSamples[i*2]);
  //ErrExit("check fuckall.csv");
  fftw_execute((*plan));
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glBegin(GL_LINE_STRIP);
  for(i = 2; i < fft_bin / 2; i++){
    complex_mag = sqrt(fft_buff[i][0]*fft_buff[i][0] +
		       fft_buff[i][1]*fft_buff[i][1]);
    glVertex3f((i-fft_bin/2.0f)/fft_bin,
	       complex_mag/(0.1f+max_value),
	       0.0);
    //printf("i: %d, cm: %f\n", i, complex_mag*1.0f);
    //glVertex3f(0.0, 0.0, 0.0);
    //glVertex3f(15, 0, 0);
    if(complex_mag > max_value){
      max_index = i;
      max_value = complex_mag;
    }
  }
  glEnd();
  glfwSwapBuffers();
  Fr = ((float)max_index * SAMPLE_RATE / pulse_size);
  distance = ((float)max_index * rr);
  printf("Frequency = %.2f Hz -----> %.2f meters\n", DecRound(Fr, 2),
	 DecRound(distance, 2));
  
  /*FILE* fout = fopen("fft_dump.csv", "w");
  for(i = 0; i < fft_bin; i++)
    fprintf(fout, "%f\n", sqrt(fft_buff[i][0]*fft_buff[i][0] +
			       fft_buff[i][1]*fft_buff[i][1]));
  ErrExit("Look at your file.");
  */



  //  ErrExit("loook at your data....");

  //  Vr = (Fd * c) / (2 * Ft); // meters/sec
  //MPH = (Vr * 3600) / (0.0254 * 12 * 5280);
  //printf("Speed = %.2f MPH. Frequency = %.2f Hz. Amplitude = %.2f\n", 
  //	 DecRound(MPH, 2), DecRound(Fd, 2), DecRound(max_value, 2));
  //printf("Speed = %.2f MPH\n", DecRound(MPH, 2));
}



void CreateScope()
{
  if(!glfwInit())
    ErrExit("Failed to initialize window.");
  
  glfwOpenWindowHint(GLFW_FSAA_SAMPLES, 4);
  glfwOpenWindowHint(GLFW_OPENGL_VERSION_MAJOR, 2);
  glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 1);
  
  if(!glfwOpenWindow(1024, 768, 0,0,0,0, 32,0, GLFW_WINDOW)){
    glfwTerminate();
    ErrExit("Failed to open GLFW window.");
  }
  
  if (glewInit() != GLEW_OK)
    ErrExit("Failed to initialize GLEW.");
  
  glfwSetWindowTitle("Radar Scope");
  glfwEnable(GLFW_STICKY_KEYS);
  glClearColor(0.0f, 0.0f, 0.4f, 0.0f);
  
  glLineWidth(2.5); 
  glColor3f(1.0, 0.0, 0.0);
}

void Artist(){
  if(!glfwInit())
    ErrExit("Failed to initialize window.");
  
  glfwOpenWindowHint(GLFW_FSAA_SAMPLES, 4);
  glfwOpenWindowHint(GLFW_OPENGL_VERSION_MAJOR, 2);
  glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 1);
  
  if(!glfwOpenWindow(1024, 768, 0,0,0,0, 32,0, GLFW_WINDOW)){
    glfwTerminate();
    ErrExit("Failed to open GLFW window.");
  }
  
  if (glewInit() != GLEW_OK)
    ErrExit("Failed to initialize GLEW.");
  
  glfwSetWindowTitle("Radar Scope");
  glfwEnable(GLFW_STICKY_KEYS);
  glClearColor(0.0f, 0.0f, 0.4f, 0.0f);
  
  glLineWidth(2.5); 
  glColor3f(1.0, 0.0, 0.0);
  glBegin(GL_LINES);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(15, 0, 0);
  glEnd();
  
  do{
    // Draw nothing, see you in tutorial 2 !
    
    // Swap buffers
    glfwSwapBuffers();
    
  } // Check if the ESC key was pressed or the window was closed
  while(glfwGetKey(GLFW_KEY_ESC) != GLFW_PRESS &&
	glfwGetWindowParam(GLFW_OPENED));
  
  // Close OpenGL window and terminate GLFW
  glfwTerminate();
  //return 0;
  
}


int main(){
  PaStreamParameters inputParameters, outputParameters;
  PaStream   *stream;
  paTestData data;
  int        c = 299792458;
  float      Fd;
  float      Ft = 2430; //Vt=2.40V
  float      Vr;
  float      MPH;

  InitAudioStreaming(&data, &inputParameters, &outputParameters);
  
  for(;;){
    switch (Prompt()){
    case 'R':
    case 'r':
      Range(&data, stream, &inputParameters);
      break;
    case 'D':
    case 'd':
      Doppler(&data, stream, &inputParameters);
      break;
    case 'W':
    case 'w':
      Record(&data, stream, &inputParameters);
      break;
    case 'P':
    case 'p':
      Play(&data, stream, &outputParameters);
      break;
    case 'A':
    case 'a':
      Artist();
      break;
    case 'S':
    case 's':
      Fd = 728.71; //Test Hz --> This would be IF
      Vr = (Fd * c) / (2 * Ft * 1000 * 1000); // meters/sec
      MPH = (Vr * 3600) / (0.0254 * 12 * 5280);

      printf("Speed = %f MPH\n", MPH);
      printf("Speed = %.4f MPH\n", DecRound(MPH, 4));
      break;
    case 'Q':
    case 'q':
      Pa_Terminate();
      printf("\n=== Your Program Works!! ===\n");
      printf("Exiting.\n\n\n");
      return c;
    }
  }
}

