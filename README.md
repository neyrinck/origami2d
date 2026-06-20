## ORIGAMI-2D

Order-ReversIng
Gravity, Apprehended
Mangling Indices

A 2D python version of ORIGAMI (Falck, Neyrinck & Szalay 2012). The ORIGAMI morphology is designed to be the number of orthogonal axes in which a particle has crossed some other particle from the initial to final conditions. In the present version, the particle x-coordinate rank (the order) is compared in rows, and the y-coordinate rank is compared in columns; if a particle is out-of-order with respect to the initial-conditions order, it is tagged. The initial ordering is assumed to be encoded in the ordering of the particle array. (In this implementation, there are 2 numpy x and y arrays, x (NxN) and y (NxN), ordered such that x increments in the 0th axis of the array and y does not change; y increments in the 1st axis of the array and x does not change.) The orderings are compared both along the Cartesian axis, and 45-degrees diagonal to it. 

<img width="1200" height="1200" alt="origami_example" src="https://github.com/user-attachments/assets/3a692bcc-b9f0-4077-9afe-9bfb03e09c4c" />

This simple comparison of initial and final orderings works most of the time, but particles can happen retain their initial ordering; here is a Lagrangian view, one pixel per particle, of the raw classification (dark purple=void; green=filament; yellow=halo). There are some apparently incorrect spots inside haloes.

<img width="640" height="480" alt="lagrangian_morph_raw" src="https://github.com/user-attachments/assets/cf1adbf5-d7cb-40d7-bbc0-2c11ba5d833c" />

With the skimage.morphology.remove_small_holes function, at least for haloes that we expect not to have holes, we can plausibly capture nearly all(?) of these underestimates of the morphology number:

<img width="640" height="480" alt="lagrangian_haloes_holesfilled" src="https://github.com/user-attachments/assets/9b96c2ec-3789-4d96-bd86-61a58125d954" />

And these can be partitioned into separate haloes with the scipy.ndimage.label function:

<img width="640" height="480" alt="group_blobs" src="https://github.com/user-attachments/assets/610e5fa2-52b3-40b3-8f38-04450c9092f2" />

Note that presently, the code does tag particles across periodic boundaries, but grouping across such boundaries does NOT happen.
