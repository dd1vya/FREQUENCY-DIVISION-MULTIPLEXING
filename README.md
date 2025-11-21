# EVALUATION OF TIME DIVISION MULTIPLEXING
### Aim:
Study of TDM pulse amplitude modulation/ demodulation with transmitter block (clock) and channel identification
information linked directly to the receivers.

### Apparatus Required:
1) Experimental kit DCL-02
2) Connecting chord
3) Power supply
4) 20 MHz dual trace oscilloscope.
   
### Theory:
In PAM, PPM the pulse is present for a short duration and for most of the time between the two pulses no signal is
present. This free space between the pulses can be occupied by pulses from other channels. This is known as Time
Division Multiplexing. Thus, time division multiplexing makes maximum utilization of the transmission channel.
Each channel to be transmitted is passed through the low pass filter. The outputs of the low pass filters are connected
to the rotating sampling switch (or) commutator.
It takes the sample from each channel per revolution and rotates at the rate of f s. Thus the sampling frequency
becomes fs the single signal composed due to multiplexing of input channels. These channels signals are then passed
through low pass reconstruction filters. If the highest signal frequency present in all the channels is fm, then by
sampling theorem, the sampling frequency fs must be such that fs≥2fm. Therefore, the time space between successive
samples from any one input will be Ts=1/fs, and Ts 1/2fm.

### Procedure:
1) Refer to the block diagram and carry out the following connections and switch setting.
2) Connect power supply in proper polarity DCL-02 and s it switch it on.
3) Connect 250hz, 500hz, 1khz and 2khz sine wave signals from the function generator to the multiplexer input
channels CH0, CH1,CH3 by means of connecting chords.
4) Connect the multiplexer output txd of the transmitted section to the de multiplexer input rxd of the receiver
section.
5) Connect the output of the receiver section CH0, CH1, CH2,CH3 to the IN0, IN1, IN2. IN3 of the filter section.
6) Connect the sampling clock TX CLK channel identification clock TXSYNC of the transmitter section to the
corresponding RXCLK, RXSYNC of the receiver section respectively.
7) Set the amplitude of the input sine wave desired.
8) Take the observations as mentioned below.
   
### Kit Diagram:
<img width="1245" height="700" alt="image" src="https://github.com/user-attachments/assets/6487a1dc-1f93-4ffc-befc-60cf907c7192" />

### Model Graph:
<img width="615" height="808" alt="image" src="https://github.com/user-attachments/assets/1d1f8d48-f08b-40fa-900c-0cb962a7e316" />

### Result:
Thus the time division multiplexing is done experimentally and output is verified.
