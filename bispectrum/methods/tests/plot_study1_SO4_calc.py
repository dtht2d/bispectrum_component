import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
folder_path = '/Users/duonghoang/Documents/GitHub/bispectrum_component/plots'  # Replace with the path to your folder

image_paths = [os.path.join(folder_path, file_name) for file_name in os.listdir(folder_path) if file_name.endswith('.png')]

# Set the number of rows and columns based on the number of images
num_images = len(image_paths)
num_rows = int(num_images ** 0.5)  # Square root of num_images rounded down
num_cols = (num_images + num_rows - 1) // num_rows  # Ceiling division
#Plot
fig, axs = plt.subplots(num_rows, num_cols, figsize=(12, 8))
for i, path in enumerate(image_paths):
    img = mpimg.imread(path)
    row = i // num_cols
    col = i % num_cols
    axs[row, col].imshow(img)
    axs[row, col].axis('off')
plt.tight_layout()
plt.show()

