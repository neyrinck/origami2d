## ORIGAMI-2D

Order-ReversIng
Gravity, Apprehended
Mangling Indices

A 2D python implementation in the same spirit as ORIGAMI (Falck, Neyrinck & Szalay 2012, MNRAS, https://ui.adsabs.harvard.edu/abs/2012ApJ...754..126F/abstract);

Order-ReversIng

Gravity, Apprehended

Mangling Indices

The ORIGAMI morphology is defined as the number of orthogonal axes in which a particle has crossed some other particle from the initial to final conditions. (This is what the original, 3D version does, https://github.com/bfalck/origami.) In the present version, a few-line, an extremely fast and simple version of that is used, which more explicitly detects "mangled indices" as in the acronym expansion. (This method misses a small fraction of crossings, as shown below, but most of these can be recovered.) This version was developed for 2D halo finding for the paper "Galaxy and Halo Root Systems: Fingerprints of Mass Assembly," https://arxiv.org/abs/2503.21015

The particle x-coordinate rank (the "order") is compared to the initial-conditions ordering in rows, and the y-coordinate rank is compared in columns; if a particle is out-of-order with respect to the initial-conditions order, it is tagged. The initial ordering is assumed to be encoded in the ordering of the particle array. (In this implementation, there are 2 numpy x and y arrays, x (NxN) and y (NxN), ordered such that x increments in the 0th axis of the array and y does not change; y increments in the 1st axis of the array and x does not change.) The orderings are compared both along the Cartesian axis, and 45-degrees diagonal to it. 

This is an example particle-tagging effort using a 2D N-body simulation, computed within the example code, with the code of Hidding (2020), https://github.com/jhidding/nbody2d

<img width="1200" height="1200" alt="origami_example" src="https://github.com/user-attachments/assets/3a692bcc-b9f0-4077-9afe-9bfb03e09c4c" />

This simple comparison of initial and final orderings tags crossed particles most of the time, but (unlike in the original ORIGAMI approach) it can fail to tag particles inside structures, because particles in a random strucure can happen to scatter back to their initial rank number in a row or column; here is a Lagrangian view, one pixel per particle, of the raw classification (dark purple=void; green=filament; yellow=halo). There are some apparently incorrect spots inside haloes.

<img width="640" height="480" alt="lagrangian_morph_raw" src="https://github.com/user-attachments/assets/cf1adbf5-d7cb-40d7-bbc0-2c11ba5d833c" />

It is less clear how to deal with this issue for filamnets, and we do not presently. But with the skimage.morphology.remove_small_holes function, we can plausibly capture nearly all(?) of these underestimates of the morphology number, by filling holes, as shown below. This involves an assumption/definition that haloes are simply-connected blobs in Lagrangian space.

<img width="640" height="480" alt="lagrangian_haloes_holesfilled" src="https://github.com/user-attachments/assets/9b96c2ec-3789-4d96-bd86-61a58125d954" />

And these can be partitioned into separate haloes with the scipy.ndimage.label function:

<img width="640" height="480" alt="group_blobs" src="https://github.com/user-attachments/assets/610e5fa2-52b3-40b3-8f38-04450c9092f2" />

Note that presently, the code does tag particles across periodic boundaries, but halo grouping and hole-filling across such boundaries does NOT happen.
