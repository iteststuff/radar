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
#define SAMPLE_TIME 0.1
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
void DopplerSample();
void Range();
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
  printf("Something went wrong: %s Exiting!", errString);
  exit(0);
}


char Prompt(){
  char c;

  printf("\nWhat would you like to do?\n");
  printf("D. Get some doppler!\n");
  printf("R. Record some audio.\n");
  printf("P. Play some audio.\n");
  printf("S. Speed Calc.\n");
  printf("W. Open window.\n");
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
    Pa_Sleep(100);
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
    Pa_Sleep(100);
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

  printf("Speed = %.2f MPH. Frequency = %.2f Hz. Amplitude = %.2f\n", 
	 DecRound(MPH, 2), DecRound(Fd, 2), DecRound(max_value, 2));
}


void Range(){
  int fstart = 2401; //Vt=2.00V
  int fstop  = 2496; //Vt=3.40V
  int bw     = fstop - fstart; //95MHz
  int tp     = 20e-3; //pulse time period
  int N      = tp * SAMPLE_RATE;
  int c      = 299792458;
  int rr     = c / (2 * bw);
  int max_range = rr * N / 2;





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
    case 'D':
    case 'd':
      Doppler(&data, stream, &inputParameters);
      break;
    case 'R':
    case 'r':
      Record(&data, stream, &inputParameters);
      break;
    case 'P':
    case 'p':
      Play(&data, stream, &outputParameters);
      break;
    case 'W':
    case 'w':
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

