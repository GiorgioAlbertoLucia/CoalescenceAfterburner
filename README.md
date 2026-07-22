Coalescence for A-nucleon nuclei
* 3 and 4 body implemented for He3, He4 (potentially Li4, but it is still not an option in the main)

How to run:
* Build with cmake:
  ```
  mkdir build && cd build
  cmake ..
  make
  ```
* Run either the He3 or He4 version (standard is 20 core multithreaded processing)
  ```
  # run He3, 100k events, radius = 1.97 fm --> it takes less then below, within few minutes
  ./he3_a3 /home/galucia/CoalescenceAfterburner/input/spectra_0_10.root ../output/he3.root 0 hProton_0_10 100000 1.97
  # run He4, 10k events, radius = 1.67 fm --> estimated time on my machine for 10k events is 353755 s ~ 100 h (CPU time). At 20 cores it is ~ 5 h
  ./li4_a4 /home/galucia/CoalescenceAfterburner/input/spectra_0_10.root ../output/he4_a4_r1p7_NEW.root 0 hProton_0_10 10000 1.67
  ```
