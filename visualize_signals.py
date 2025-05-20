import matplotlib.pyplot as plt
from scikit_build_example import _core

# Parametry sygnału
freq = 440           # Hz
sample_rate = 44100  # Hz
samples = 1000       # liczba próbek

# Generowanie sygnałów
sine = _core.generate_sine_wave(freq, sample_rate, samples)
cosine = _core.generate_cosine_wave(freq, sample_rate, samples)
square = _core.generate_square_wave(freq, sample_rate, samples)
sawtooth = _core.generate_sawtooth_wave(freq, sample_rate, samples)

# Oś czasu w sekundach
time = [i / sample_rate for i in range(samples)]

# Tworzenie wykresów
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(time, sine)
plt.title("Sine Wave")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.subplot(2, 2, 2)
plt.plot(time, cosine)
plt.title("Cosine Wave")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.subplot(2, 2, 3)
plt.plot(time, square)
plt.title("Square Wave")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.subplot(2, 2, 4)
plt.plot(time, sawtooth)
plt.title("Sawtooth Wave")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.tight_layout()
plt.show()
