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

filtered_sine = _core.filter_1d(sine, [0.2, 0.6, 0.2])

# Przykładowy 2D sygnał – np. szum lub gradient
image = [[(i + j) % 255 for j in range(100)] for i in range(100)]

# Przykładowy filtr dolnoprzepustowy (rozmycie)
kernel = [
    [1/9, 1/9, 1/9],
    [1/9, 1/9, 1/9],
    [1/9, 1/9, 1/9]
]

filtered_image = _core.filter_2d(image, kernel)

dft_result = _core.dft(sine)
reconstructed = _core.idft(dft_result)

# Tworzenie wykresów
plt.figure(figsize=(12, 8))

plt.subplot(4, 2, 1)
plt.plot(time, sine)
plt.title("Sine Wave")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.subplot(4, 2, 2)
plt.plot(time, cosine)
plt.title("Cosine Wave")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.subplot(4, 2, 3)
plt.plot(time, square)
plt.title("Square Wave")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.subplot(4, 2, 4)
plt.plot(time, sawtooth)
plt.title("Sawtooth Wave")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.subplot(4, 2, 5)
plt.plot(filtered_sine)
plt.title("Filtrowana fala sinusoidalna")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.subplot(4, 2, 6)
plt.imshow(filtered_image, cmap='gray')
plt.title("Po filtracji")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.subplot(4, 2, 7)
plt.plot([abs(x) for x in dft_result])
plt.title("DFT - amplituda (moduł)")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.subplot(4, 2, 8)
plt.plot(reconstructed)
plt.title("Zrekonstruowana fala po IDFT")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.tight_layout()
plt.show()
