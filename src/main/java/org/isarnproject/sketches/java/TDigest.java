/*
Copyright 2016 Erik Erlandson

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

package org.isarnproject.sketches.java;

public class TDigest {
    public double compression;
    public double[] centers;
    public double[] masses;
    public double[] cumulative;

    public TDigest() {
        int n = 20;
        centers = new double[n];
        masses = new double[n];
        cumulative = new double[1 + n];
        for (int j = 0; j < centers.length; ++j) {
            centers[j] = 0.0;
            masses[j] = 0.0;
            cumulative[j] = 0.0;
        }
        cumulative[n] = 0.0;
    }

    public static final int lsb(int j) {
        return j & (-j);
    }

    public double sumFT(int j) {
        j += 1;
        double s = 0.0;
        while (j > 0) {
            s += cumulative[j];
            j -= lsb(j);
        }
        return s;
    }

    public void addFT(int j, double x) {
        int n = cumulative.length;
        j += 1;
        while (j < n) {
            cumulative[j] += x;
            j += lsb(j);
        }
    }
}
