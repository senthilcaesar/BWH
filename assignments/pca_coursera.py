import os
import cv2
import numpy as np
import matplotlib.pyplot as plt


def load_images_from_directory(directory_path):
    images = []
    for filename in os.listdir(directory_path):
        print(filename)
        file_path = os.path.join(directory_path, filename)
        img = cv2.imread(file_path, cv2.IMREAD_GRAYSCALE)
        images.append(img)
    return images

def center_data(Y):
    mean_vector = np.mean(Y, axis=0)
    mean_matrix = Y - mean_vector    
    return mean_matrix

def get_cov_matrix(X):
    cov_matrix = np.dot(X.T, X)
    cov_matrix = cov_matrix/(X.shape[0]-1)    
    return cov_matrix

def perform_PCA(X, eigenvecs, k):
    V = eigenvecs[:,:k]
    Xred = np.dot(X, V)
    return Xred

def reconstruct_image(Xred, eigenvecs):
    X_reconstructed = Xred.dot(eigenvecs[:,:Xred.shape[1]].T)
    return X_reconstructed

directory_path = '/Users/sq566/BWH/assignments/CatDog'
imgs = load_images_from_directory(directory_path)

height, width = imgs[0].shape
print(f'\nYour dataset has {len(imgs)} images of size {height}x{width} pixels\n')


plt.imshow(imgs[15], cmap='gray')

imgs_flatten = np.array([im.reshape(-1) for im in imgs])
print(f'imgs_flatten shape: {imgs_flatten.shape}')

X = center_data(imgs_flatten)
plt.imshow(X[15].reshape(64,64), cmap='gray')

cov_matrix = get_cov_matrix(X)
print(f'Covariance matrix shape: {cov_matrix.shape}')


eigenvals, eigenvecs = np.linalg.eig(cov_matrix)
eigenvals = eigenvals.real
eigenvecs = eigenvecs.real
idx = np.argsort(np.abs(eigenvals))[::-1]
sorted_eigenvalues = eigenvals[idx]
sorted_eigenvectors = eigenvecs[:, idx]

# Select the top 55 components
eigenvals = sorted_eigenvalues[:55]
eigenvecs = sorted_eigenvectors[:, :55]
print(f'Ten largest eigenvalues: \n{eigenvals[:10]}')

fig, ax = plt.subplots(4,4, figsize=(20,20))
for n in range(4):
    for k in range(4):
        ax[n,k].imshow(eigenvecs[:,n*4+k].reshape(height,width), cmap='gray')
        ax[n,k].set_title(f'component number {n*4+k+1}')
        
        
        
        
Xred2 = perform_PCA(X, eigenvecs,2)
print(f'Xred2 shape: {Xred2.shape}')     
        


x = Xred2[:, 0]
y = Xred2[:, 1]

# Create the plot
plt.figure(figsize=(10, 6))
plt.scatter(x, y, c='blue', label='Data points', alpha=0.5)
for i, (x_val, y_val) in enumerate(zip(x, y)):
    plt.text(x_val, y_val, str(i), color='black', fontsize=8)
plt.grid(False)
plt.show()


fig, ax = plt.subplots(1,3, figsize=(15,5))
ax[0].imshow(imgs[14], cmap='gray')
ax[0].set_title('Image 14')
ax[1].imshow(imgs[25], cmap='gray')
ax[1].set_title('Image 25')
ax[2].imshow(imgs[43], cmap='gray')
ax[2].set_title('Image 43')
plt.suptitle('Similar cats')



fig, ax = plt.subplots(1,3, figsize=(15,5))
ax[0].imshow(imgs[14], cmap='gray')
ax[0].set_title('Image 14')
ax[1].imshow(imgs[45], cmap='gray')
ax[1].set_title('Image 45')
ax[2].imshow(imgs[49], cmap='gray')
ax[2].set_title('Image 49')
plt.suptitle('Different cats')



Xred1 = perform_PCA(X, eigenvecs,1) # reduce dimensions to 1 component
Xred5 = perform_PCA(X, eigenvecs, 5) # reduce dimensions to 5 components
Xred10 = perform_PCA(X, eigenvecs, 10) # reduce dimensions to 10 components
Xred20 = perform_PCA(X, eigenvecs, 20) # reduce dimensions to 20 components
Xred30 = perform_PCA(X, eigenvecs, 30) # reduce dimensions to 30 components
Xrec1 = reconstruct_image(Xred1, eigenvecs) # reconstruct image from 1 component
Xrec5 = reconstruct_image(Xred5, eigenvecs) # reconstruct image from 5 components
Xrec10 = reconstruct_image(Xred10, eigenvecs) # reconstruct image from 10 components
Xrec20 = reconstruct_image(Xred20, eigenvecs) # reconstruct image from 20 components
Xrec30 = reconstruct_image(Xred30, eigenvecs) # reconstruct image from 30 components



fig, ax = plt.subplots(2,3, figsize=(22,15))
ax[0,0].imshow(imgs[14], cmap='gray')
ax[0,0].set_title('original', size=20)
ax[0,1].imshow(Xrec1[14].reshape(height,width), cmap='gray')
ax[0,1].set_title('reconstructed from 1 components', size=20)
ax[0,2].imshow(Xrec5[14].reshape(height,width), cmap='gray')
ax[0,2].set_title('reconstructed from 5 components', size=20)
ax[1,0].imshow(Xrec10[14].reshape(height,width), cmap='gray')
ax[1,0].set_title('reconstructed from 10 components', size=20)
ax[1,1].imshow(Xrec20[14].reshape(height,width), cmap='gray')
ax[1,1].set_title('reconstructed from 20 components', size=20)
ax[1,2].imshow(Xrec30[14].reshape(height,width), cmap='gray')
ax[1,2].set_title('reconstructed from 30 components', size=20)

fig = plt.figure()
explained_variance = 100 * eigenvals/sum(eigenvals)
plt.plot(np.arange(1,56), explained_variance)
plt.xlabel("Component Number")
plt.ylabel("\u03BB")

fig = plt.figure()
explained_cum_variance = np.cumsum(explained_variance/100)
plt.plot(np.arange(1,56), explained_cum_variance)
plt.xlabel("Component Number")
plt.axhline(y=0.95, color='r')



Xred35 = perform_PCA(X, eigenvecs, 35) # reduce dimensions to 35 components
Xrec35 = reconstruct_image(Xred35, eigenvecs) # reconstruct image from 35 components

fig, ax = plt.subplots(4,2, figsize=(15,28))
ax[0,0].imshow(imgs[0], cmap='gray')
ax[0,0].set_title('original', size=20)
ax[0,1].imshow(Xrec35[0].reshape(height, width), cmap='gray')
ax[0,1].set_title('Reconstructed', size=20)

ax[1,0].imshow(imgs[15], cmap='gray')
ax[1,0].set_title('original', size=20)
ax[1,1].imshow(Xrec35[15].reshape(height, width), cmap='gray')
ax[1,1].set_title('Reconstructed', size=20)

ax[2,0].imshow(imgs[32], cmap='gray')
ax[2,0].set_title('original', size=20)
ax[2,1].imshow(Xrec35[32].reshape(height, width), cmap='gray')
ax[2,1].set_title('Reconstructed', size=20)

ax[3,0].imshow(imgs[54], cmap='gray')
ax[3,0].set_title('original', size=20)
ax[3,1].imshow(Xrec35[54].reshape(height, width), cmap='gray')
ax[3,1].set_title('Reconstructed', size=20)

















































