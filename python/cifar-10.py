import tensorflow as tf
from tensorflow.keras import layers, models
import matplotlib.pyplot as plt
import pandas as pd
import os

# Create results directory if it doesn't exist
os.makedirs('results', exist_ok=True)

# 1. Load and prepare CIFAR-10
(x_train, y_train), (x_test, y_test) = tf.keras.datasets.cifar10.load_data()
x_train, x_test = x_train / 255.0, x_test / 255.0

# 2. Define the model
model = models.Sequential([
    layers.Conv2D(32, (3, 3), activation='relu', input_shape=(32, 32, 3)),
    layers.MaxPooling2D((2, 2)),
    layers.Conv2D(64, (3, 3), activation='relu'),
    layers.MaxPooling2D((2, 2)),
    layers.Conv2D(64, (3, 3), activation='relu'),
    layers.Flatten(),
    layers.Dense(64, activation='relu'),
    layers.Dense(10)
])

# 3. Compile
model.compile(optimizer='adam',
              loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
              metrics=['accuracy'])

# 4. Train and save the history
print("Training started... This will take a few minutes.")
history = model.fit(x_train, y_train, epochs=10, 
                    validation_data=(x_test, y_test))

# 5. EVALUATE AND CREATE THE TABLE (CSV)
print("\nGenerating results table...")
# Create a dictionary of the training metrics
results_data = {
    'Epoch': range(1, 11),
    'Accuracy': history.history['accuracy'],
    'Loss': history.history['loss'],
    'Val_Accuracy': history.history['val_accuracy'],
    'Val_Loss': history.history['val_loss']
}
# Convert to a Pandas DataFrame and save as CSV
df = pd.DataFrame(results_data)
df.to_csv('results/training_results.csv', index=False)
print(" Table saved to: results/training_results.csv")

# 6. PLOT THE FIGURES
print("\nGenerating plots...")
plt.figure(figsize=(12, 4))

# Plot Accuracy
plt.subplot(1, 2, 1)
plt.plot(results_data['Epoch'], results_data['Accuracy'], label='Training Accuracy')
plt.plot(results_data['Epoch'], results_data['Val_Accuracy'], label='Validation Accuracy')
plt.title('Training and Validation Accuracy')
plt.xlabel('Epochs')
plt.ylabel('Accuracy')
plt.legend()
plt.grid(True)

# Plot Loss
plt.subplot(1, 2, 2)
plt.plot(results_data['Epoch'], results_data['Loss'], label='Training Loss')
plt.plot(results_data['Epoch'], results_data['Val_Loss'], label='Validation Loss')
plt.title('Training and Validation Loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.grid(True)

# Save the plot
plt.tight_layout()
plt.savefig('results/cifar_training_plot.png', dpi=300)
plt.show()
print(" Plot saved to: results/cifar_training_plot.png")

# 7. Final Test Accuracy
test_loss, test_acc = model.evaluate(x_test, y_test, verbose=2)
print(f'\nFinal Test accuracy: {test_acc}')
