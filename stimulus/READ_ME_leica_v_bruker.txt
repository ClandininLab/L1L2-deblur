There are two Stimulus files in 2p-stim-code, and the correct one needs to be named "Stimulus" based on which setup you're using.

By default, there should be a file titled "Stimulus" which should be for the Leica, and another file titled "Stimulus_bruker" which should be for the Bruker. If this is the case and you are imaging on the Leica, you can use 2p-stim-code as is.

If you would like to image on the Bruker, you need to:
1. change "Stimulus" to "Stimulus_leica"
2. change "Stimulus_bruker" to "Stimulus"
3. image
4. revert to the original names after you are done!

(If, instead, you see a file upon startup titled "Stimulus" and another file titled "Stimulus_leica", then someone has forgotten to rename the files after they finished imaging on the Bruker. In order to image on the Leica, you need to:
1. change "Stimulus" to "Stimulus_bruker"
2. change "Stimulus_leica" to "Stimulus".)
