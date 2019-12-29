c = get_config()
c.ExecutePreprocessor.timeout = 300
c.ExecutePreprocessor.kernel_name='python3'
c.NbConvertApp.notebooks = ["../atmosphere/Atmosphere.ipynb",
                            "../atmosphere/AtmospherePhaseScreenGradient.ipynb",
                            "../LMMSE/LaserTomography.ipynb",
                            "../shackHartmann/Shack-HartmannWavefrontSensor.ipynb",
                            "../shackHartmann/Shack-HartmannNoisePropagation.ipynb",
                            "../shackHartmann/LtaoShWfs.ipynb",
                            "../rayTracing/GMT-Ray-Tracing.ipynb",
                            "../segmentPistonSensor/Segment-Piston-Sensor.ipynb"]
